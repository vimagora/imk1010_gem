"""Tests for strain construction.

Several of these tests assert a *disagreement* between the model and the
published phenotype. That is deliberate. The disagreement is the main finding
of Part A, so it should fail loudly if a future change to the code or the model
makes it quietly go away.
"""

import pytest

from imk1010_gem import analysis, strains
from imk1010_gem.model_io import load_model


@pytest.fixture(scope="module")
def model():
    return load_model()


@pytest.fixture(scope="module")
def built(model):
    """IMK1010 and its gene mapping, built once."""
    return strains.build_imk1010(model)


@pytest.fixture(scope="module")
def repaired(model):
    substrates = list(analysis.DEVALK_SUBSTRATES.values())
    return strains.build_repaired_imk1010(model, substrates)


# --- the published gene list ------------------------------------------------


def test_gene_list_matches_the_published_counts():
    """de Valk et al. report 23 transporter and 9 hydrolase genes."""
    rows = strains.read_deleted_genes()
    categories = [row["category"] for row in rows]
    assert categories.count("transporter") == 23
    assert categories.count("hydrolase") == 9
    assert len(rows) == 32


def test_three_genes_are_absent_from_yeast9(built):
    """S288C, which Yeast9 is built on, lacks the CEN.PK MAL2 locus."""
    _, mapping = built
    assert len(mapping.mapped) == 29
    assert set(mapping.unmapped) == {"HXT12", "MAL21", "MAL22"}


# --- direction handling ------------------------------------------------------


def test_export_only_channel_is_not_counted_as_uptake(model):
    """r_1172 is written cytosol -> extracellular and is irreversible."""
    channel = model.reactions.get_by_id(strains.GLYCEROL_CHANNEL)
    assert channel.lower_bound == 0
    assert not strains.can_import(channel, "glycerol")


def test_symporter_is_counted_as_uptake(model):
    symporter = model.reactions.get_by_id("r_1171")
    assert strains.can_import(symporter, "glycerol")


# --- the finding -------------------------------------------------------------


def test_sucrose_and_maltose_uptake_are_genuinely_blocked(model, built):
    _, mapping = built
    for substrate in ["sucrose", "maltose"]:
        assert strains.surviving_uptake_genes(model, substrate, mapping.orfs) == {}


def test_hexose_uptake_survives_the_published_deletions(model, built):
    """HXT17 and CIN10 are in Yeast9's glucose rule but were never deleted."""
    _, mapping = built
    survivors = strains.surviving_uptake_genes(model, "D-glucose", mapping.orfs)
    assert "r_1166" in survivors
    assert set(survivors["r_1166"]) == {"YNR072W", "YDR387C"}  # HXT17, CIN10


def test_faithful_imk1010_wrongly_grows_on_glucose(built):
    """The Part A disagreement: published phenotype is no growth."""
    strain, _ = built
    assert analysis.grows(analysis.growth_on(strain, "D-glucose"))


def test_repaired_strain_cannot_grow_on_glucose(repaired):
    strain, _, extra = repaired
    assert not analysis.grows(analysis.growth_on(strain, "D-glucose"))
    assert "YNR072W" in extra["D-glucose"]  # HXT17 had to be deleted


def test_repaired_strain_predicts_no_glycerol_growth(repaired):
    """Blocking the hexose leak also blocks glycerol.\n\n    This matches Nummela et al. (2026), whose empty-vector control did not\n    grow on 1% glycerol, but contradicts de Valk et al. (2022), who saw\n    growth on 2% glycerol after 8 days. The two models bracket the two\n    observations; see scripts/run_part_a.py.\n    """
    strain, _, _ = repaired
    assert not analysis.grows(analysis.growth_on(strain, "glycerol"))


def test_bidirectional_channel_restores_glycerol_without_unblocking_hexoses(repaired):
    """The resolution: Fps1 is a channel, so it should conduct both ways."""
    strain, _, _ = repaired
    with strain:
        strains.allow_glycerol_channel_uptake(strain)
        assert analysis.grows(analysis.growth_on(strain, "glycerol"))
        assert not analysis.grows(analysis.growth_on(strain, "D-glucose"))
        assert not analysis.grows(analysis.growth_on(strain, "D-fructose"))


def test_ethanol_growth_is_predicted_by_every_variant(model, built, repaired):
    """Ethanol needs no sugar transporter, so all variants should agree."""
    faithful, _ = built
    repaired_strain, _, _ = repaired
    for strain in [model, faithful, repaired_strain]:
        assert analysis.grows(analysis.growth_on(strain, "ethanol"))
