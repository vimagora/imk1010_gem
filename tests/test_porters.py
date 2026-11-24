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
