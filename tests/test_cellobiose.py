"""Tests for the cellobiose platform strain (Part C)."""

import pytest

from imk1010_gem import analysis, cellobiose, strains
from imk1010_gem.model_io import load_model

SUBSTRATES = ["D-glucose", "D-fructose", "D-mannose", "D-glucitol", "xylitol", "glycerol"]


@pytest.fixture(scope="module")
def model():
    return load_model()


@pytest.fixture(scope="module")
def repaired(model):
    strain, _, _ = strains.build_repaired_imk1010(model, SUBSTRATES)
    return strain


def test_cellobiose_is_absent_from_yeast9(model):
    """Part C is an extension, so the starting point must be empty."""
    assert not [m for m in model.metabolites if "cellobiose" in m.name.lower()]


def test_hydrolysis_is_mass_balanced(repaired):
    """C12H22O11 + H2O -> 2 C6H12O6. If this does not balance, nothing downstream means anything."""
    with repaired as strain:
        _, inside = cellobiose.add_cellobiose_pools(strain)
        reaction = cellobiose.add_intracellular_bgl(strain, inside=inside)
        assert reaction.check_mass_balance() == {}


def test_hydrolysis_involves_no_redox_cofactor(repaired):
    """A beta-glucosidase is a hydrolase: no NAD(P)H, so no cofactor imbalance."""
    with repaired as strain:
        _, inside = cellobiose.add_cellobiose_pools(strain)
        reaction = cellobiose.add_intracellular_bgl(strain, inside=inside)
        names = {m.name.upper() for m in reaction.metabolites}
        assert not names & {"NAD", "NADH", "NADP(+)", "NADPH", "NADP", "FAD", "FADH2"}


def test_platform_strain_cannot_grow_without_a_transporter(model):
    """The point of the design: the enzyme alone is not enough."""
    strain = cellobiose.build_imk1010_ncbgl(model, SUBSTRATES)
    assert not analysis.grows(analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME))


def test_adding_a_transporter_restores_growth(model):
    strain = cellobiose.build_imk1010_ncbgl(model, SUBSTRATES)
    with strain:
        cellobiose.add_cellobiose_transporter(strain, "CtA")
        assert analysis.grows(analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME))


@pytest.mark.parametrize("proton_symport", [True, False])
def test_growth_does_not_depend_on_the_transport_mechanism(model, proton_symport):
    """Mechanism is unconfirmed for CtA; it must not change the conclusion."""
    strain = cellobiose.build_imk1010_ncbgl(model, SUBSTRATES)
    with strain:
        cellobiose.add_cellobiose_transporter(strain, "CtA", proton_symport=proton_symport)
        assert analysis.grows(analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME))


def test_proton_symport_costs_growth_but_only_a_little(model):
    """One symport event carries two hexoses, so the per-hexose cost is halved."""
    rates = {}
    for symport in [True, False]:
        strain = cellobiose.build_imk1010_ncbgl(model, SUBSTRATES)
        with strain:
            cellobiose.add_cellobiose_transporter(strain, "CtA", proton_symport=symport)
            rates[symport] = analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME)
    assert rates[True] < rates[False]
    assert (1 - rates[True] / rates[False]) < 0.10


def test_secreted_bgl_rescues_wild_type_but_not_imk1010(model, repaired):
    """Localisation is the whole design: glucose released outside cannot get in."""
    with model as wild_type:
        cellobiose.add_cellobiose_pools(wild_type)
        cellobiose.add_secreted_bgl(wild_type)
        assert analysis.grows(analysis.growth_on(wild_type, cellobiose.CELLOBIOSE_NAME))

    with repaired as strain:
        cellobiose.add_cellobiose_pools(strain)
        cellobiose.add_secreted_bgl(strain)
        assert not analysis.grows(analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME))


def test_constructs_data_file_records_sources():
    """Every construct must carry an accession and a citation it can be traced to."""
    rows = cellobiose.read_constructs()
    constructs = {r["construct"] for r in rows}
    assert {"NcBGL (GH1-1)", "CtA", "secretedBGL"} <= constructs

    ncbgl = next(r for r in rows if r["construct"].startswith("NcBGL"))
    assert ncbgl["localisation"] == "intracellular"
    assert ncbgl["accession"] == "XP_011393758.1"
    assert ncbgl["gene_or_locus"] == "NCU00130"

    cta = next(r for r in rows if r["construct"] == "CtA")
    assert cta["accession"] == "MH648002"
    assert cta["gene_or_locus"] == "An12g09270"
    # The name clash with A. nidulans CltA/CltB is a real trap; keep it recorded.
    assert "CltA" in cta["notes"]

    # Every construct taken from the literature cites where it was characterised.
    for row in rows:
        if row["source_organism"].startswith("("):
            continue  # control constructs have no source paper
        assert row["characterised_in"].strip()
        assert row["doi_or_pmid"].strip()
