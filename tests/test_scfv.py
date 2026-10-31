"""scfv: crystal chain positions to scFv positions and back, on the trastuzumab scFv (VH 120, linker 15, VL 107).

The offset is tested against a second, independent construction: edit the VH and the VL strings on their own and join
them with the linker, which does not use the layout at all.
"""

import pytest

from design_pipeline import scfv

PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    "GGGGSGGGGSGGGGS"
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
VH, VL = PARENT.split(scfv.LINKER)
HEAVY_CRYSTAL = VH + "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC"
LIGHT_CRYSTAL = VL + "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"


@pytest.fixture
def layout():
    return scfv.layout_from_parent(PARENT)


def test_the_parent_splits_into_the_two_domains_and_the_linker(layout):
    assert (len(VH), len(VL), len(PARENT)) == (120, 107, 242)
    assert layout.offsets == {"B": 0, "A": 135} and layout.lengths == {"B": 120, "A": 107}
    assert layout.linker == (120, 15)
    assert layout.domain("B") == VH and layout.domain("A") == VL


def test_the_heavy_chain_keeps_its_numbers_and_the_light_chain_shifts_by_135(layout):
    assert [scfv.scfv_position(layout, "B", p) for p in (1, 97, 120)] == [1, 97, 120]
    assert [scfv.scfv_position(layout, "A", p) for p in (1, 27, 107)] == [136, 162, 242]


@pytest.mark.parametrize("chain, position", [("B", 0), ("B", 121), ("B", 220), ("A", 108), ("A", 214), ("A", 0), ("A", -3)])
def test_a_constant_domain_or_missing_position_has_no_counterpart(layout, chain, position):
    with pytest.raises(ValueError, match="outside the variable domain"):
        scfv.scfv_position(layout, chain, position)


def test_an_unknown_chain_is_refused(layout):
    with pytest.raises(ValueError, match="not in the layout"):
        scfv.scfv_position(layout, "H", 5)


def test_every_variable_position_round_trips_and_the_linker_has_no_chain_position(layout):
    for chain, length in (("B", 120), ("A", 107)):
        for position in range(1, length + 1):
            assert scfv.chain_position(layout, scfv.scfv_position(layout, chain, position)) == (chain, position)
    for index in range(121, 136):
        with pytest.raises(ValueError, match="in the linker"):
            scfv.chain_position(layout, index)
    for index in (0, 243):
        with pytest.raises(ValueError, match="outside the sequence"):
            scfv.chain_position(layout, index)


def test_a_mutation_in_each_domain_lands_where_editing_the_domains_alone_puts_it(layout):
    """The independent construction: change VH and VL as strings, then rejoin. No layout involved."""
    heavy, light = (97, VH[96], "A"), (27, VL[26], "E")
    expected_vh = VH[:96] + "A" + VH[97:]
    expected_vl = VL[:26] + "E" + VL[27:]
    result = scfv.apply_mutations(layout, [("B", heavy[0], heavy[1], heavy[2]), ("A", light[0], light[1], light[2])])
    assert result == expected_vh + scfv.LINKER + expected_vl
    assert len(result) == 242 and len(scfv.substitutions(PARENT, result)) == 2


def test_apply_id_reads_the_variant_id_format_and_ignores_the_order_of_the_tokens(layout):
    a = scfv.apply_id(layout, f"B:{VH[96]}97A+A:{VL[26]}27E")
    b = scfv.apply_id(layout, f"A:{VL[26]}27E+B:{VH[96]}97A")
    assert a == b and a != PARENT


def test_the_wild_type_is_checked_against_the_parent_so_a_shift_is_an_error(layout):
    """Light-chain position 27 given the residue that sits at scFv position 27 of the heavy chain, as if there were no offset."""
    wrong = PARENT[26]
    assert wrong != VL[26]
    with pytest.raises(ValueError, match="does not match the parent"):
        scfv.apply_mutations(layout, [("A", 27, wrong, "E")])


@pytest.mark.parametrize(
    "mutations, message",
    [
        ([("B", 121, "A", "G")], "outside the variable domain"),
        ([("B", 5, VH[4], "A"), ("B", 5, VH[4], "G")], "two mutations at B:5"),
        ([("B", 5, VH[4], VH[4])], "not a substitution"),
        ([("B", 5, VH[4], "X")], "not a substitution"),
        ([("B", 5, VH[4], "AC")], "not a substitution"),
    ],
)
def test_a_candidate_that_cannot_be_applied_is_refused_with_the_reason(layout, mutations, message):
    with pytest.raises(ValueError, match=message):
        scfv.apply_mutations(layout, mutations)


def test_a_variant_goes_back_to_the_crystal_numbering(layout):
    mutations = [("B", 97, VH[96], "A"), ("B", 103, VH[102], "Y"), ("A", 27, VL[26], "E"), ("A", 107, VL[106], "R"), ("B", 1, VH[0], "D")]
    variant = scfv.apply_mutations(layout, mutations)
    assert sorted(scfv.crystal_mutations(layout, variant)) == sorted(mutations)
    assert len(scfv.substitutions(PARENT, variant)) == 5


def test_a_change_in_the_linker_is_refused_on_the_way_back(layout):
    variant = PARENT[:125] + ("A" if PARENT[125] != "A" else "G") + PARENT[126:]
    with pytest.raises(ValueError, match="in the linker"):
        scfv.crystal_mutations(layout, variant)


def test_substitutions_only_so_a_different_length_raises():
    with pytest.raises(ValueError, match="only substitutions"):
        scfv.substitutions(PARENT, PARENT + "A")
    assert scfv.substitutions(PARENT, PARENT) == []


def test_a_light_chain_first_scfv_puts_the_offset_on_the_heavy_chain():
    flipped = VL + scfv.LINKER + VH
    layout = scfv.layout_from_parent(flipped, first="A", second="B")
    assert layout.offsets == {"A": 0, "B": 122} and scfv.scfv_position(layout, "A", 1) == 1
    assert scfv.scfv_position(layout, "B", 1) == 123 and scfv.scfv_position(layout, "B", 120) == 242


@pytest.mark.parametrize(
    "bad, message",
    [
        (VH + VL, "exactly once"),
        (VH + scfv.LINKER + VL + scfv.LINKER + VL, "exactly once"),
        (scfv.LINKER + VL, "one side of it is empty"),
        (VH + scfv.LINKER, "one side of it is empty"),
        (PARENT[:10] + "Z" + PARENT[11:], "not amino acids"),
    ],
)
def test_a_parent_that_does_not_have_the_shape_is_refused(bad, message):
    with pytest.raises(ValueError, match=message):
        scfv.layout_from_parent(bad)


def test_two_domains_need_two_chains():
    with pytest.raises(ValueError, match="two different chains"):
        scfv.layout_from_parent(PARENT, first="B", second="B")


def test_verify_parent_accepts_a_crystal_whose_variable_domains_are_the_scfvs(layout, monkeypatch):
    monkeypatch.setattr(scfv, "chain_sequence", lambda structure, chain: {"B": HEAVY_CRYSTAL, "A": LIGHT_CRYSTAL}[chain])
    scfv.verify_parent(layout, "1N8Z.cif")


def test_verify_parent_names_the_chain_and_position_where_the_crystal_differs(layout, monkeypatch):
    changed = LIGHT_CRYSTAL[:49] + ("G" if LIGHT_CRYSTAL[49] != "G" else "A") + LIGHT_CRYSTAL[50:]
    monkeypatch.setattr(scfv, "chain_sequence", lambda structure, chain: {"B": HEAVY_CRYSTAL, "A": changed}[chain])
    with pytest.raises(ValueError, match=r"chain A of 1N8Z.cif .*first difference at position 50"):
        scfv.verify_parent(layout, "1N8Z.cif")


def test_verify_parent_refuses_a_crystal_chain_shorter_than_the_domain(layout, monkeypatch):
    monkeypatch.setattr(scfv, "chain_sequence", lambda structure, chain: {"B": VH[:100], "A": LIGHT_CRYSTAL}[chain])
    with pytest.raises(ValueError, match=r"chain B .*first difference at position 101 \(crystal \(end of chain\)"):
        scfv.verify_parent(layout, "x.cif")
