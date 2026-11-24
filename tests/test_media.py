"""Tests for medium setup.

These guard the one assumption every later result rests on: that when we say
"growth on glycerol", the model really had glycerol and nothing else as a
carbon source.
"""

import pytest

from imk1010_gem import media
from imk1010_gem.model_io import load_model


@pytest.fixture(scope="module")
def model():
    """Load Yeast9 once and share it; reading the 12 MB SBML is slow."""
    return load_model()


def test_model_loads_with_expected_size(model):
    # Guards against silently analysing a different model than we documented.
    assert len(model.reactions) == 4105
    assert len(model.metabolites) == 2748
    assert len(model.genes) == 1143


def test_objective_is_growth(model):
    objective = [r for r in model.reactions if r.objective_coefficient != 0]
    assert len(objective) == 1
    assert objective[0].id == "r_2111"
    assert objective[0].name == "growth"


def test_shipped_medium_has_glucose_as_only_carbon_source(model):
    assert media.open_carbon_sources(model) == ["D-glucose"]


def test_wild_type_grows_on_glucose(model):
    with model:
        media.set_carbon_source(model, "D-glucose")
        assert model.optimize().objective_value > 1e-6


@pytest.mark.parametrize(
    "substrate", ["D-glucose", "D-fructose", "D-mannose", "glycerol", "xylitol", "D-glucitol"]
)
def test_exactly_one_carbon_source_is_open(model, substrate):
    """The core invariant of media.set_carbon_source."""
    with model:
        exchange = media.set_carbon_source(model, substrate)
        open_now = media.open_carbon_sources(model)
        assert len(open_now) == 1, f"expected 1 carbon source, got {open_now}"
        assert exchange.lower_bound < 0


def test_medium_change_is_reverted_by_context_manager(model):
    with model:
        media.set_carbon_source(model, "glycerol")
        assert media.open_carbon_sources(model) == ["glycerol"]
    assert media.open_carbon_sources(model) == ["D-glucose"]


def test_unknown_substrate_raises_rather_than_guessing(model):
    with pytest.raises(KeyError):
        media.find_exchange(model, "D-sorbitol")  # stored as 'D-glucitol' in Yeast9


def test_per_carbon_basis_scales_by_formula(model):
    """Six carbon atoms of glycerol (C3) is 2 mmol of glycerol."""
    with model:
        exchange = media.set_carbon_source(model, "glycerol", uptake=6.0, per_carbon=True)
        assert exchange.lower_bound == pytest.approx(-2.0)
