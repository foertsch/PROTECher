"""Step 1: obtaining the complex, offline."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from design_pipeline import Complex, RunFailedError, get_structure
from design_pipeline.adapters._pdb import FIXTURE_CHAINS, write_ca_pdb
from design_pipeline.structures import (
    chain_ids,
    chain_sequence,
    load,
    read_fasta,
)


@pytest.fixture
def structure(tmp_path: Path) -> Path:
    return write_ca_pdb(tmp_path / "source" / "1abc.pdb", FIXTURE_CHAINS)


@pytest.fixture
def cif_structure(tmp_path: Path, structure: Path) -> Path:
    """The same complex as mmCIF: what the co-folding models write, and the
    only format a very large assembly is published in."""
    from biotite.structure.io import pdbx

    out = pdbx.CIFFile()
    pdbx.set_structure(out, load(structure))
    path = tmp_path / "source" / "1abc.cif"
    out.write(str(path))
    return path


def test_the_download_path_runs_from_a_local_source(tmp_path, structure) -> None:
    """`source` replaces the HTTP GET, which is how this runs with no network."""
    complex_ = get_structure(
        tmp_path / "runs", source=structure, chains=["H", "L", "A"]
    )

    assert isinstance(complex_, Complex)
    assert complex_.path.is_file()
    assert complex_.path != structure  # copied into the run directory
    assert chain_ids(complex_.path) == ("H", "L", "A")
    assert complex_.config["model_name"] == "pdb_download"


def test_a_missing_chain_is_caught_before_anything_is_scored(tmp_path, structure):
    """The labelling trap: 1N8Z calls its heavy chain B. Finding out now beats
    finding out after 4,000 variants have been scored at the wrong positions."""
    with pytest.raises(RunFailedError, match="does not have chain"):
        get_structure(tmp_path / "runs", source=structure, chains=["H", "Z"])


def test_a_complex_is_frozen_and_its_representative_is_the_first_path(tmp_path):
    complex_ = get_structure(tmp_path / "runs", model_name="dummy_complex")

    assert len(complex_) == 1
    assert complex_.path is complex_.paths[0]
    with pytest.raises(Exception):
        complex_.paths = ()  # type: ignore[misc]


def test_the_run_is_cached_and_resumes(tmp_path, structure) -> None:
    first = get_structure(tmp_path / "runs", source=structure)
    second = get_structure(tmp_path / "runs", source=structure)

    assert first.path == second.path
    assert first.config["cache_key"] == second.config["cache_key"]
    # One run directory, not two: the key is the resolved config plus the
    # content of the input.
    assert len([p for p in (tmp_path / "runs").iterdir() if p.is_dir()]) == 1


def test_sequences_are_written_to_a_fasta_for_a_predictor(tmp_path) -> None:
    """Every predictor in step 1 takes the same input: one FASTA record with
    the chains joined by ':'. Checked here with the fixture folder, because the
    real ones need a GPU."""
    get_structure(
        tmp_path / "runs",
        model_name="dummy_structures",
        antibody_seq="EVQLVE:DIQMTQ",
        antigen_seq="TQVCTG",
    )
    fasta = tmp_path / "runs" / "inputs" / "complex.fasta"

    assert fasta.read_text().splitlines()[1] == "EVQLVE:DIQMTQ:TQVCTG"


def test_a_structure_with_no_accession_and_no_source_is_an_error(tmp_path) -> None:
    with pytest.raises(ValueError, match="pdb_id"):
        get_structure(tmp_path / "runs")


def test_chain_sequence_indexes_the_chain_as_given(tmp_path, structure) -> None:
    """`position` in the variant table is 1-based into this string, not the
    PDB residue number."""
    assert chain_sequence(structure, "H") == FIXTURE_CHAINS["H"]
    assert chain_sequence(structure, "A") == FIXTURE_CHAINS["A"]


def test_a_cif_source_keeps_its_own_format(tmp_path, cif_structure) -> None:
    """`file_format` picks what to ask RCSB for; it does not rename a source.

    Naming the copy .pdb was silent here and fatal later: `load` dispatches on
    the suffix, so the mmCIF came back out through the PDB reader at the first
    step that read a sequence off it.
    """
    complex_ = get_structure(tmp_path / "runs", source=cif_structure)

    assert complex_.path.suffix == ".cif"
    assert chain_ids(complex_.path) == ("H", "L", "A")


def test_a_compressed_source_is_refused_rather_than_relabelled(tmp_path, structure):
    """The RCSB download links everyone pastes are gzipped."""
    archive = tmp_path / "source" / "1abc.pdb.gz"
    archive.write_bytes(gzip.compress(structure.read_bytes()))

    with pytest.raises(ValueError, match="Decompress"):
        get_structure(tmp_path / "runs", source=archive)


def test_file_format_names_a_source_that_does_not_name_itself(tmp_path, structure):
    """A suffix that says nothing about the format is the one case where
    `file_format` still decides."""
    plain = tmp_path / "source" / "mailed_to_me.txt"
    plain.write_bytes(structure.read_bytes())

    complex_ = get_structure(tmp_path / "runs", source=plain, file_format="pdb")

    assert complex_.path.suffix == ".pdb"
    assert chain_ids(complex_.path) == ("H", "L", "A")


def test_load_refuses_a_file_it_cannot_place(tmp_path, structure) -> None:
    unknown = tmp_path / "1abc.structure"
    unknown.write_bytes(structure.read_bytes())

    with pytest.raises(ValueError, match="cannot tell the format"):
        load(unknown)



def test_a_fasta_on_disk_replaces_the_sequence_arguments(tmp_path) -> None:
    """The sequences you already have in a file, rather than pasted into the
    call. One record per chain is the layout a database or an alignment tool
    hands you."""
    source = tmp_path / "my_complex.fasta"
    source.write_text(">heavy\nEVQLVE\n>light\nDIQMTQ\n>antigen\nTQVCTG\n")

    get_structure(tmp_path / "runs", model_name="dummy_structures", fasta=source)

    normalised = tmp_path / "runs" / "inputs" / "complex.fasta"
    assert normalised.read_text().splitlines()[1] == "EVQLVE:DIQMTQ:TQVCTG"


def test_a_fasta_and_the_sequence_arguments_agree_on_one_input(tmp_path) -> None:
    """Either layout normalises to the same input file, so a campaign started
    from a FASTA and one started from the arguments share a cache key."""
    source = tmp_path / "joined.fasta"
    source.write_text(">complex\nEVQLVE:DIQMTQ:TQVCTG\n")

    from_file = get_structure(
        tmp_path / "a", model_name="dummy_structures", fasta=source
    )
    from_args = get_structure(
        tmp_path / "b",
        model_name="dummy_structures",
        antibody_seq="EVQLVE:DIQMTQ",
        antigen_seq="TQVCTG",
    )

    assert from_file.config["cache_key"] == from_args.config["cache_key"]


def test_a_wrapped_multi_record_fasta_reads_as_its_chains(tmp_path) -> None:
    """A FASTA body is wrapped at 60 or 80 columns; the wrap is not a chain
    break."""
    source = tmp_path / "wrapped.fasta"
    source.write_text(">heavy\nEVQL\nVESG\n\n>light desc here\nDIQM\n")

    assert read_fasta(source) == ["EVQLVESG", "DIQM"]


def test_sequences_and_a_fasta_together_are_refused(tmp_path) -> None:
    source = tmp_path / "x.fasta"
    source.write_text(">complex\nEVQLVE\n")

    with pytest.raises(ValueError, match="not both"):
        get_structure(
            tmp_path / "runs",
            model_name="dummy_structures",
            fasta=source,
            antibody_seq="EVQLVE",
        )


def test_sequences_handed_to_the_download_adapter_say_so(tmp_path) -> None:
    """`pdb_download` is the default model_name, so this is the easy mistake:
    the error used to complain about a missing pdb_id instead."""
    source = tmp_path / "x.fasta"
    source.write_text(">complex\nEVQLVE\n")

    with pytest.raises(ValueError, match="takes no sequences"):
        get_structure(tmp_path / "runs", fasta=source)


def test_a_structure_handed_to_fasta_is_caught(tmp_path, structure) -> None:
    with pytest.raises(ValueError, match="holds no FASTA records"):
        get_structure(
            tmp_path / "runs", model_name="dummy_structures", fasta=structure
        )


def test_a_gapped_alignment_is_not_a_sequence(tmp_path) -> None:
    source = tmp_path / "aligned.fasta"
    source.write_text(">heavy\nEVQL--GG\n")

    with pytest.raises(ValueError, match="not a residue"):
        read_fasta(source)
