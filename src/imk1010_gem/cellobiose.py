"""Part C: a cellobiose-utilising platform strain, IMK1010 + NcBGL.

The logic of the strain
-----------------------
Saccharomyces cerevisiae cannot use cellobiose: it has neither a cellobiose
transporter nor a beta-glucosidase. Two things must therefore be added, and
*which* of them is limiting is the experiment.

IMK1010 + NcBGL (an intracellular beta-glucosidase, and nothing else) should
**not** grow on cellobiose, because the disaccharide cannot cross the
membrane. Adding a cellobiose transporter should restore growth. That makes
the strain a growth-based reporter for cellobiose *transport*, which is how
such transporters are characterised experimentally, and it mirrors Part B:
a dead platform strain, then one heterologous protein at a time.

Why the beta-glucosidase must be intracellular
----------------------------------------------
A *secreted* beta-glucosidase hydrolyses cellobiose outside the cell and
releases glucose into the medium. Wild-type yeast grows on that happily.
IMK1010 cannot: it has no hexose transporter left, so the glucose stays
outside. This is modelled explicitly (``add_secreted_bgl``) because it is a
sharp demonstration of what the Sugar0 background is for, and a reminder that
"beta-glucosidase" alone is an under-specified construct -- several
well-known BGL1 enzymes are secreted.

What this part does not do
--------------------------
It does not rank transporters. Yeast9 has no affinities, so every cellobiose
transporter with the same stoichiometry yields an identical prediction; the
17 cellobiose transporters compiled in Xu et al. (2024) Table S1 would give 17
identical rows. The question a genome-scale model can answer is whether
transport alone suffices, or whether a pathway is missing too.
"""

from __future__ import annotations

import csv
import pathlib

import cobra

from .strains import build_repaired_imk1010

CONSTRUCTS_CSV = (
    pathlib.Path(__file__).resolve().parents[2] / "data" / "cellobiose_constructs.csv"
)

# Cellobiose is absent from Yeast9 entirely, so both pools are created here.
# beta-1,4-linked glucose dimer: C12H22O11.
CELLOBIOSE_FORMULA = "C12H22O11"
CELLOBIOSE_NAME = "cellobiose"


def read_constructs(path: pathlib.Path | str = CONSTRUCTS_CSV) -> list[dict[str, str]]:
    """Read the heterologous constructs and their sources."""
    with pathlib.Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _metabolite(model: cobra.Model, name: str, compartment: str) -> cobra.Metabolite:
    """Look a metabolite up by name and compartment, failing loudly."""
    matches = [
        metabolite
        for metabolite in model.metabolites
        if metabolite.name.strip().lower() == name.lower()
        and metabolite.compartment == compartment
    ]
    if len(matches) != 1:
        raise KeyError(
            f"expected exactly one {name!r} in compartment {compartment!r}, "
            f"found {[m.id for m in matches]}"
        )
    return matches[0]


def add_cellobiose_pools(model: cobra.Model) -> tuple[cobra.Metabolite, cobra.Metabolite]:
    """Create extracellular and cytosolic cellobiose, plus an exchange reaction.

    THIS EXTENDS YEAST9. Cellobiose is not in the published model, so this is
    an addition of ours and is reported as such, never as a model result.

    Both pools carry the same name; only the extracellular one gets an
    exchange, so medium lookups by name remain unambiguous.
    """
    outside = cobra.Metabolite(
        "s_cellobiose_e", name=CELLOBIOSE_NAME, formula=CELLOBIOSE_FORMULA,
        charge=0, compartment="e",
    )
    inside = cobra.Metabolite(
        "s_cellobiose_c", name=CELLOBIOSE_NAME, formula=CELLOBIOSE_FORMULA,
        charge=0, compartment="c",
    )
    exchange = cobra.Reaction("EX_cellobiose_e", name="cellobiose exchange")
    exchange.add_metabolites({outside: -1})
    exchange.bounds = (0.0, 1000.0)  # closed until a medium opens it
    model.add_reactions([exchange])
    return outside, inside


def add_intracellular_bgl(
    model: cobra.Model, gene: str = "NcBGL", inside: cobra.Metabolite | None = None
) -> cobra.Reaction:
    """Add the intracellular beta-glucosidase: cellobiose + H2O -> 2 D-glucose.

    This is a hydrolysis, not an oxidation. No NAD(P)H is involved and there is
    no redox cost or cofactor imbalance -- unlike, say, a cellobiose
    dehydrogenase, or the xylose reductase/xylitol dehydrogenase pair.

    Modelled after N. crassa GH1-1 (NCU00130, XP_011393758.1), which is
    intracellular. Kinetic constants are not needed: FBA uses stoichiometry,
    not rates.
    """
    inside = inside or _metabolite(model, CELLOBIOSE_NAME, "c")
    glucose_c = _metabolite(model, "D-glucose", "c")
    water_c = _metabolite(model, "H2O", "c")

    reaction = cobra.Reaction(
        f"{gene}_cellobiose_hydrolysis",
        name=f"cellobiose hydrolysis by {gene} (intracellular beta-glucosidase)",
    )
    reaction.add_metabolites({inside: -1, water_c: -1, glucose_c: 2})
    reaction.bounds = (0.0, 1000.0)
    reaction.gene_reaction_rule = gene
    model.add_reactions([reaction])
    return reaction


def add_secreted_bgl(model: cobra.Model, gene: str = "secretedBGL") -> cobra.Reaction:
    """The negative-design control: hydrolysis *outside* the cell.

    Included to show why localisation is the whole point. In wild type this
    restores growth; in IMK1010 it cannot, because the glucose released into
    the medium has no transporter to enter by.
    """
    outside = _metabolite(model, CELLOBIOSE_NAME, "e")
    glucose_e = _metabolite(model, "D-glucose", "e")
    water_e = _metabolite(model, "H2O", "e")

    reaction = cobra.Reaction(
        f"{gene}_cellobiose_hydrolysis",
        name=f"cellobiose hydrolysis by {gene} (secreted beta-glucosidase)",
    )
    reaction.add_metabolites({outside: -1, water_e: -1, glucose_e: 2})
    reaction.bounds = (0.0, 1000.0)
    reaction.gene_reaction_rule = gene
    model.add_reactions([reaction])
    return reaction


def add_cellobiose_transporter(
    model: cobra.Model, gene: str = "CtA", proton_symport: bool = True
) -> cobra.Reaction:
    """Add a cellobiose transporter.

    Args:
        gene: transporter name, used as the heterologous gene ID.
        proton_symport: whether uptake is coupled to a proton. The mechanism
            is NOT confirmed for CtA, so this is an assumption and both
            settings are reported in scripts/run_part_c.py. It does not change
            growth/no-growth; it changes the ATP cost, and therefore yield.

    A detail worth noticing: one symport event delivers *two* hexose units, so
    the proton cost per hexose is half that of proton-coupled glucose uptake.
    That is the published energetic argument for cellobiose co-fermentation.
    """
    outside = _metabolite(model, CELLOBIOSE_NAME, "e")
    inside = _metabolite(model, CELLOBIOSE_NAME, "c")

    reaction = cobra.Reaction(
        f"t_{gene}_cellobiose",
        name=f"cellobiose transport via {gene}"
        + (" (H+ symport)" if proton_symport else " (facilitated)"),
    )
    stoichiometry: dict[cobra.Metabolite, float] = {outside: -1, inside: 1}
    if proton_symport:
        stoichiometry[_metabolite(model, "H+", "e")] = -1
        stoichiometry[_metabolite(model, "H+", "c")] = 1
    reaction.add_metabolites(stoichiometry)
    reaction.bounds = (0.0, 1000.0)
    reaction.gene_reaction_rule = gene
    model.add_reactions([reaction])
    return reaction


def build_imk1010_ncbgl(
    model: cobra.Model, substrates: list[str], gene: str = "NcBGL"
) -> cobra.Model:
    """The platform strain: IMK1010 (repaired) + cellobiose pools + NcBGL.

    No transporter. It should not grow on cellobiose; that is the point.
    """
    strain, _, _ = build_repaired_imk1010(model, substrates)
    strain.id = f"IMK1010_{gene}"
    _, inside = add_cellobiose_pools(strain)
    add_intracellular_bgl(strain, gene=gene, inside=inside)
    return strain
