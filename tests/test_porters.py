"""Tests for adding the Glp porters."""

import pytest

from imk1010_gem import analysis, porters, strains
from imk1010_gem.model_io import load_model


@pytest.fixture(scope="module")
def model():
    return load_model()


@pytest.fixture(scope="module")
def background(model):
    """IMK1010 repaired: the background the porters were expressed in."""
    substrates = ["D-glucose", "D-fructose", "D-mannose", "D-glucitol", "xylitol", "glycerol"]
    strain, _, _ = strains.build_repaired_imk1010(model, substrates)
    return strain


def test_published_specificity_is_read_from_the_data_file():
    specificity = porters.read_glp_specificity()
    assert set(specificity) == {"GlpA", "GlpB", "GlpC"}
    # All three carry the polyols; only GlpC carries hexoses.
    for porter in ["GlpA", "GlpB", "GlpC"]:
        assert "glycerol" in specificity[porter]
    assert "D-glucose" in specificity["GlpC"]
    assert "D-glucose" not in specificity["GlpA"]


def test_background_strain_grows_on_nothing(background):
    """Matches Nummela's empty-vector control; without this there is no signal."""
    for substrate in ["D-glucose", "D-fructose", "D-mannose", "D-glucitol",
                      "xylitol", "glycerol"]:
        assert not analysis.grows(analysis.growth_on(background, substrate)), substrate


def test_adding_a_porter_restores_growth_on_its_substrates(background):
    specificity = porters.read_glp_specificity()
    with background as strain:
        porters.add_porter(strain, "GlpA", specificity["GlpA"])
        for substrate in ["glycerol", "D-glucitol", "xylitol"]:
            assert analysis.grows(analysis.growth_on(strain, substrate)), substrate


def test_porter_does_not_restore_substrates_it_does_not_carry(background):
    """The discriminating test: GlpA/GlpB were not reported to carry hexoses."""
    specificity = porters.read_glp_specificity()
    with background as strain:
        porters.add_porter(strain, "GlpB", specificity["GlpB"])
        for hexose in ["D-glucose", "D-fructose", "D-mannose"]:
            assert not analysis.grows(analysis.growth_on(strain, hexose)), hexose


def test_glpc_restores_hexoses_too(background):
    specificity = porters.read_glp_specificity()
    with background as strain:
        porters.add_porter(strain, "GlpC", specificity["GlpC"])
        for hexose in ["D-glucose", "D-fructose", "D-mannose"]:
            assert analysis.grows(analysis.growth_on(strain, hexose)), hexose


def test_porter_reactions_are_reverted_by_the_context_manager(background):
    before = len(background.reactions)
    with background as strain:
        porters.add_porter(strain, "GlpA", ["glycerol"])
        assert len(strain.reactions) == before + 1
    assert len(background.reactions) == before


def test_mannitol_is_not_representable_without_an_explicit_extension():
    assert porters.unrepresentable_substrates(["D-mannitol"]) == ["D-mannitol"]
    assert analysis.YEAST9_NAMES["D-mannitol"] == ""


def test_mannitol_gap_is_transport_not_pathway(background):
    """Yeast9 has mannitol -> D-fructose (r_4183) but no way to take it up.

    Adding only the uptake side restores growth, so the gap is transport.
    """
    extended = background.copy()
    porters.add_extracellular_mannitol(extended)
    with extended as strain:
        porters.add_mannitol_porter(strain, "GlpA")
        assert analysis.grows(analysis.growth_on(strain, "Mannitol"))


# --- the published negatives, which are what make Part B informative ---------


def test_glpa_and_glpb_hexose_results_are_published_negatives():
    """Tested and did not grow -- not 'unreported'.

    Encoding these as absent results would hide the discriminating half of
    Part B: a model that said yes to every added transporter would score all
    the positives and none of these.
    """
    import csv

    with open(porters.GLP_CSV, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for porter in ["GlpA", "GlpB"]:
        for sugar in ["D-glucose", "D-fructose", "D-mannose"]:
            match = [r for r in rows if r["porter"] == porter and r["substrate"] == sugar]
            assert len(match) == 1, f"{porter}/{sugar}"
            assert match[0]["transported"] == "no"


def test_galactose_was_carried_by_no_porter():
    import csv

    with open(porters.GLP_CSV, newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle) if r["substrate"] == "D-galactose"]
    assert len(rows) == 3
    assert all(r["transported"] == "no" for r in rows)


def test_ethanol_is_the_published_positive_control():
    """Every strain grew on ethanol, which is what shows the strains were viable."""
    import csv

    with open(porters.GLP_CSV, newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle) if r["substrate"] == "ethanol"]
    assert len(rows) == 3
    assert all(r["transported"] == "yes" for r in rows)


def test_panel_contains_both_positives_and_negatives():
    """A panel of only positives could be matched by a trivially permissive model."""
    import csv

    with open(porters.GLP_CSV, newline="", encoding="utf-8") as handle:
        calls = [r["transported"] for r in csv.DictReader(handle)]
    assert calls.count("yes") == 18
    assert calls.count("no") == 9
    assert "not reported" not in calls
