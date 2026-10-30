"""Target selection and the per-target protocol, offline.

The network functions are thin wrappers; what is tested is the part that
decides what gets benchmarked - the validation, the SAbDab filters, the
deduplication and the cutoff rules - because a quiet mistake in any of them
turns the whole exercise into a memorisation test.
"""

from __future__ import annotations

import io
import math
import warnings
from pathlib import Path

import pytest

from design_pipeline import benchmark, metrics
from design_pipeline.benchmark import BenchmarkError, Target

SUMMARY_COLUMNS = (
    "pdb\tHchain\tLchain\tantigen_chain\tantigen_type\tantigen_name\tdate\t"
    "resolution\tmethod\n"
)
ROWS = (
    ("8aaa", "H", "L", "C", "protein", "spike", "2023-05-01", "2.10", "X-RAY DIFFRACTION"),
    ("8bbb", "H", "L", "C", "protein", "spike", "2023-04-01", "2.50", "X-RAY DIFFRACTION"),
    ("8ccc", "H", "L", "C", "protein", "lysozyme", "2023-03-01", "1.90", "X-RAY DIFFRACTION"),
    ("8ddd", "H", "", "C", "protein", "nanobody target", "2023-02-01", "2.00", "X-RAY DIFFRACTION"),
    ("8eee", "H", "L", "C", "peptide", "a peptide", "2023-01-01", "2.00", "X-RAY DIFFRACTION"),
    ("8fff", "H", "L", "H", "protein", "itself", "2023-01-01", "2.00", "X-RAY DIFFRACTION"),
    ("8ggg", "H", "L", "C", "protein", "too old", "2019-01-01", "2.00", "X-RAY DIFFRACTION"),
    ("8hhh", "H", "L", "C", "protein", "cryo", "2023-06-01", "3.80", "ELECTRON MICROSCOPY"),
)


@pytest.fixture
def summary():
    import pandas as pd

    text = SUMMARY_COLUMNS + "\n".join("\t".join(row) for row in ROWS)
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

GOOD_ENTRY = {
    "rcsb_entry_info": {
        "resolution_combined": [2.1],
        "polymer_entity_count_protein": 3,
    },
    "rcsb_accession_info": {"initial_release_date": "2003-01-14T00:00:00Z"},
}
GOOD_SEQUENCES = {
    "B": "EVQLVESGGGLVQPGGSLRLSCAAS",
    "A": "DIQMTQSPSSLSASVGDRVTITCRAS",
    "C": "TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQ",
}
TARGET = Target(pdb_id="1N8Z", heavy="B", light="A", antigen="C", tier="smoke")


def test_a_good_entry_passes() -> None:
    assert benchmark.check_target(TARGET, GOOD_ENTRY, GOOD_SEQUENCES) == []


def test_a_low_resolution_entry_is_rejected() -> None:
    entry = {"rcsb_entry_info": {**GOOD_ENTRY["rcsb_entry_info"], "resolution_combined": [3.4]}}
    failures = benchmark.check_target(TARGET, entry, GOOD_SEQUENCES)

    assert any("resolution" in reason for reason in failures)


def test_an_entry_with_no_resolution_is_rejected() -> None:
    """No resolution means it is not a crystal structure, whatever else it is."""
    entry = {"rcsb_entry_info": {"polymer_entity_count_protein": 3}}

    assert any("resolution" in r for r in benchmark.check_target(TARGET, entry, GOOD_SEQUENCES))


def test_too_few_protein_entities_is_rejected() -> None:
    """Both variable domains and a protein antigen have to be resolved; two
    entities is a nanobody or an unbound Fab."""
    entry = {
        "rcsb_entry_info": {"resolution_combined": [2.0], "polymer_entity_count_protein": 2}
    }
    failures = benchmark.check_target(TARGET, entry, GOOD_SEQUENCES)

    assert any("polymer entities" in reason for reason in failures)


def test_a_chain_the_entry_does_not_have_is_rejected() -> None:
    """The hardcoded smoke-tier chain assignments are not trusted; this is the
    check that catches one that is wrong."""
    failures = benchmark.check_target(TARGET, GOOD_ENTRY, {"X": "ACDEF" * 10})

    assert any("not in the entry" in reason for reason in failures)


def test_an_antigen_that_is_really_an_antibody_chain_is_rejected() -> None:
    sequences = {**GOOD_SEQUENCES, "C": GOOD_SEQUENCES["B"]}
    failures = benchmark.check_target(TARGET, GOOD_ENTRY, sequences)

    assert any("not the antigen" in reason for reason in failures)


def test_a_peptide_antigen_is_rejected() -> None:
    sequences = {**GOOD_SEQUENCES, "C": "ACDEFGHIK"}
    failures = benchmark.check_target(TARGET, GOOD_ENTRY, sequences)

    assert any("peptide" in reason for reason in failures)


def test_validation_drops_rather_than_raises() -> None:
    """A failed assertion prints and drops the entry; the run continues on
    whatever survived."""
    kept, dropped = benchmark.validate_smoke_tier(
        [TARGET, Target(pdb_id="9ZZZ", heavy="X", light="Y", antigen="Z")],
        fetch=lambda pdb_id: GOOD_ENTRY,
        sequences=lambda pdb_id: GOOD_SEQUENCES,
    )

    assert [t.pdb_id for t in kept] == ["1N8Z"]
    assert kept[0].resolution == 2.1
    assert kept[0].date == "2003-01-14"
    assert [t.pdb_id for t, _ in dropped] == ["9ZZZ"]


def test_an_unreachable_entry_is_a_drop_not_a_crash() -> None:
    def fetch(pdb_id):
        raise BenchmarkError("could not reach RCSB")

    kept, dropped = benchmark.validate_smoke_tier([TARGET], fetch=fetch, sequences=None)

    assert kept == []
    assert "could not reach RCSB" in dropped[0][1][0]


def test_the_smoke_tier_is_the_three_the_spec_names() -> None:
    assert [t.pdb_id for t in benchmark.SMOKE_TIER] == ["1VFB", "1N8Z", "3HFM"]
    assert all(t.tier == "smoke" for t in benchmark.SMOKE_TIER)
    # 1N8Z labels trastuzumab's heavy chain B, which is the trap.
    her2 = next(t for t in benchmark.SMOKE_TIER if t.pdb_id == "1N8Z")
    assert her2.heavy == "B" and her2.light == "A"


# --------------------------------------------------------------------------
# RCSB parsing
# --------------------------------------------------------------------------


def test_the_entry_fasta_is_parsed_per_chain() -> None:
    text = (
        ">1N8Z_1|Chains A, B|Immunoglobulin|Mus musculus\nQVQLQQ\n"
        ">1N8Z_2|Chain C|HER2|Homo sapiens\nTQVCTG\n"
    )
    assert benchmark.parse_rcsb_fasta(text) == {
        "A": "QVQLQQ",
        "B": "QVQLQQ",
        "C": "TQVCTG",
    }


def test_an_auth_renamed_chain_keeps_the_author_id() -> None:
    """`A[auth H]` is how RCSB writes a renamed asym id, and the author id is
    the one a PDB file - and so every metric here - actually uses."""
    parsed = benchmark.parse_rcsb_fasta(">x_1|Chain A[auth H]|Ig|Mus\nQVQL\n")

    assert parsed == {"H": "QVQL"}


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2023-05-01", "2023-05-01"),
        ("05/01/23", "2023-05-01"),
        ("05/01/2023", "2023-05-01"),
        ("", None),
        ("not a date", None),
    ],
)
def test_sabdab_dates_are_read_in_either_layout(text, expected) -> None:
    """SAbDab has written this column both ways over the years."""
    parsed = benchmark.parse_date(text)

    assert (parsed.isoformat() if parsed else None) == expected


# --------------------------------------------------------------------------
# Training cutoffs
# --------------------------------------------------------------------------


def test_alphafold3_has_a_recorded_cutoff() -> None:
    cutoff, source, verified = benchmark.training_cutoff("alphafold3")

    assert cutoff == "2021-09-30"
    assert verified is True
    assert "paper" in source.lower()


def test_an_undocumented_predictor_refuses_to_invent_one() -> None:
    """A date guessed here would be indistinguishable from one that was
    checked, and getting it wrong puts training data into the honest tier."""
    cutoff, source, verified = benchmark.training_cutoff("boltz2")

    assert cutoff is None
    assert verified is False
    assert "model card" in source


def test_the_release_minus_six_months_fallback() -> None:
    assert benchmark.cutoff_from_release("2024-06-01") == "2023-12-01"
    with pytest.raises(ValueError, match="not a date"):
        benchmark.cutoff_from_release("whenever")


# --------------------------------------------------------------------------
# SAbDab selection
# --------------------------------------------------------------------------


def test_the_filters_run_in_order_and_report_each_count(summary) -> None:
    targets, log = benchmark.select_sabdab(summary, "2021-09-30", limit=6)

    assert log[0] == "8 rows in the summary file"
    assert "7 with a protein antigen" in log[1]
    assert "6 with a heavy chain" in log[2]
    assert "5 whose antigen chain is not an antibody chain" in log[3]
    assert "4 solved by X-ray" in log[4]
    assert "3 deposited after 2021-09-30" in log[5]
    assert [t.pdb_id for t in targets] == ["8AAA", "8CCC"]


def test_the_most_recent_survive_and_carry_their_metadata(summary) -> None:
    targets, _ = benchmark.select_sabdab(summary, "2021-09-30", limit=6)

    assert targets[0].pdb_id == "8AAA"
    assert targets[0].date == "2023-05-01"
    assert targets[0].resolution == 2.10
    assert targets[0].tier == "honest"
    assert targets[0].chains == ("H", "L", "C")


def test_the_limit_truncates(summary) -> None:
    targets, log = benchmark.select_sabdab(summary, "2021-09-30", limit=1)

    assert len(targets) == 1
    assert log[-1] == "1 selected, most recent first"


def test_a_stricter_cutoff_only_shrinks_the_pool(summary) -> None:
    few, _ = benchmark.select_sabdab(summary, "2023-04-15", limit=6)
    many, _ = benchmark.select_sabdab(summary, "2021-09-30", limit=6)

    assert {t.pdb_id for t in few} <= {t.pdb_id for t in many}


def test_a_bad_cutoff_is_an_error(summary) -> None:
    with pytest.raises(ValueError, match="not a date"):
        benchmark.select_sabdab(summary, "soon", limit=6)


def test_antigens_are_deduplicated_by_sequence_when_sequences_are_available(summary):
    """Without this the set fills with near-identical SARS-CoV-2 RBD complexes
    and six entries are one data point wearing six hats."""
    spike = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHS"
    lysozyme = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGS"
    shared = {
        "8aaa": {"C": spike},
        "8bbb": {"C": spike[:-2] + "QQ"},  # a near-identical spike
        "8ccc": {"C": lysozyme},
    }
    targets, log = benchmark.select_sabdab(
        summary, "2021-09-30", limit=6, sequences=lambda pdb: shared.get(pdb, {})
    )

    assert any("identity" in line for line in log)
    # The two spikes collapse to one; the lysozyme is its own cluster.
    assert len(targets) == 2
    assert {t.pdb_id for t in targets} == {"8AAA", "8CCC"}


def test_without_sequences_it_falls_back_to_the_antigen_name(summary) -> None:
    targets, log = benchmark.select_sabdab(summary, "2021-09-30", limit=6)

    assert any("antigen_name" in line and "coarser" in line for line in log)
    assert len(targets) == 2


def test_identity_clustering_is_single_linkage() -> None:
    clusters = benchmark.cluster_by_identity(
        {
            "a": "ACDEFGHIKLMNPQRSTVWY",
            "b": "ACDEFGHIKLMNPQRSTVWW",  # one substitution from a
            "c": "WWWWWWWWWWWWWWWWWWWW",
        },
        threshold=0.70,
    )

    assert sorted(len(cluster) for cluster in clusters) == [1, 2]
    assert any({"a", "b"} == set(cluster) for cluster in clusters)


def test_a_summary_file_without_the_expected_columns_is_refused(tmp_path) -> None:
    path = tmp_path / "wrong.tsv"
    path.write_text("a\tb\n1\t2\n")

    with pytest.raises(BenchmarkError, match="missing column"):
        benchmark.read_sabdab(path)


def test_a_local_summary_copy_is_read(tmp_path, summary) -> None:
    path = tmp_path / "summary.tsv"
    path.write_text(SUMMARY_COLUMNS + "\n".join("\t".join(row) for row in ROWS))

    assert len(benchmark.read_sabdab(path)) == len(ROWS)


# --------------------------------------------------------------------------
# The per-target protocol
# --------------------------------------------------------------------------


def test_the_fixture_targets_look_like_deposited_entries(tmp_path) -> None:
    made = benchmark.fixture_targets(tmp_path, 3)

    assert len(made) == 3
    for target, source in made:
        assert target.tier == "fixture"
        assert target.chains == ("B", "A", "C")  # relabelled, as a real entry is
        assert source.is_file()
    # A different antigen per target, so each one gets its own ipTM and poses.
    from design_pipeline.structures import chain_sequence

    lengths = {len(chain_sequence(source, "C")) for _, source in made}
    assert len(lengths) == 3


def test_the_reference_is_relabelled_and_its_sequences_read_back(tmp_path) -> None:
    """The prediction is made from the residues that were actually resolved, so
    a disordered loop cannot turn into a length mismatch at the comparison."""
    (target, source), *_ = benchmark.fixture_targets(tmp_path / "fix", 1)
    normalized, (heavy, light, antigen) = benchmark.reference_structure(
        target, output_dir=tmp_path / "runs", reference_source=source
    )

    from design_pipeline.structures import chain_ids

    assert chain_ids(normalized) == ("H", "L", "A")
    assert len(heavy) > 20 and len(light) > 20 and len(antigen) > 20


def test_a_target_is_predicted_scored_and_recorded(tmp_path) -> None:
    (target, source), *_ = benchmark.fixture_targets(tmp_path / "fix", 1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row = benchmark.evaluate_target(
            target,
            output_dir=tmp_path / "runs",
            model_name="dummy_predictor",
            num_seeds=4,
            reference_source=source,
        )

    assert row["pdb_id"] == target.pdb_id
    assert row["tier"] == "fixture"
    assert row["n_seeds"] == 4
    assert 0.0 <= row["epitope_recall"] <= 1.0
    assert row["interface_rmsd"] == row["interface_rmsd"]  # not NaN
    assert row["confidence_iptm"] is not None
    # The same flat schema every other step writes.
    assert metrics.load(tmp_path / "runs" / target.pdb_id)["epitope_recall"] == (
        row["epitope_recall"]
    )


def test_a_missing_dockq_is_noted_rather_than_fatal(tmp_path) -> None:
    (target, source), *_ = benchmark.fixture_targets(tmp_path / "fix", 1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row = benchmark.evaluate_target(
            target,
            output_dir=tmp_path / "runs",
            model_name="dummy_predictor",
            num_seeds=2,
            reference_source=source,
        )

    try:
        import DockQ  # noqa: F401
    except ModuleNotFoundError:
        assert math.isnan(row["dockq"])
        assert "metrics" in row["dockq_note"]
    else:  # pragma: no cover - only where the extra is installed
        assert row["dockq"] == row["dockq"]
