"""Adding the Aspergillus niger Glp transporters to the model.

Nummela et al. (2026) expressed GlpA, GlpB and GlpC in IMK1010 and read out
growth on solid medium after six days. All three restored growth on glycerol,
D-mannitol, D-sorbitol and D-xylitol; GlpC additionally restored growth on
D-fructose, D-glucose and D-mannose. The empty-vector control grew on none of
them.

So the model question is qualitative and matches the experiment's readout: with
a porter's published substrates made transportable, does IMK1010 grow?

What we assume, and why it is the weakest point of Part B
---------------------------------------------------------
Each porter is added as facilitated transport, ``substrate[e] <=> substrate[c]``,
with no proton coupling. Nummela et al. do not report the transport mechanism,
so this is our assumption, chosen because:

* the Glp family are aquaglyceroporin-like channels/facilitators, and
* proton symport would add an ATP cost per molecule, which would change
  predicted *yields*; with no measured yields to compare against, assuming the
  cheaper mechanism avoids inventing a cost.

It is stated here rather than buried, because for growth/no-growth it happens
not to matter (either mechanism permits uptake), but for any future flux or
yield comparison it would matter a great deal.

A second limitation: FBA has no affinities. Nummela needed 8% glucose and 8%
fructose to see GlpC-dependent growth but only 1% for the polyols, which
suggests low affinity for hexoses. The model cannot represent that difference
at all -- in FBA a transporter either exists or it does not.
"""

from __future__ import annotations

import csv
import pathlib

import cobra

from . import analysis

GLP_CSV = pathlib.Path(__file__).resolve().parents[2] / "data" / "nummela2026_glp.csv"

# Mannitol exists in Yeast9 only as a cytosolic metabolite (s_3805), with no
# extracellular form and no exchange reaction, so the model cannot take it up
# under any circumstances. The catabolic step is present (r_4183,
# D-Mannitol:NAD+ 2-oxidoreductase, mannitol + NAD <=> D-fructose + NADH),
# so only the uptake side is missing.
CYTOSOLIC_MANNITOL = "s_3805"


def read_glp_specificity(path: pathlib.Path | str = GLP_CSV) -> dict[str, list[str]]:
    """Map each porter to the substrates it was published to transport."""
    specificity: dict[str, list[str]] = {}
    with pathlib.Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["transported"].strip().lower() == "yes":
                specificity.setdefault(row["porter"], []).append(row["substrate"])
    return specificity


def add_extracellular_mannitol(model: cobra.Model) -> tuple[cobra.Metabolite, cobra.Reaction]:
    """Give Yeast9 an extracellular mannitol pool and an exchange reaction.

    THIS EXTENDS THE MODEL BEYOND WHAT YEAST9 SHIPS. It is applied only to a
    separate, clearly-labelled strain variant, and it is reported as an
    extension rather than a result.

    The point of the extension is to separate two different explanations for
    "Yeast9 cannot grow on mannitol": a missing transporter, or a missing
    pathway. Since the catabolic step is already present, adding the uptake
    side tests whether the pathway alone suffices.

    Formula and charge are copied from the cytosolic metabolite rather than
    typed in, so the two pools cannot disagree.
    """
    cytosolic = model.metabolites.get_by_id(CYTOSOLIC_MANNITOL)

    extracellular = cobra.Metabolite(
        "s_mannitol_e",
        name=cytosolic.name,
        formula=cytosolic.formula,
        charge=cytosolic.charge,
        compartment="e",
    )
    exchange = cobra.Reaction("EX_mannitol_e", name="mannitol exchange")
    exchange.add_metabolites({extracellular: -1})
    exchange.bounds = (0.0, 1000.0)  # closed for uptake until a medium opens it
    model.add_reactions([exchange])
    return extracellular, exchange


def add_porter(
    model: cobra.Model,
    porter: str,
    substrates: list[str],
    reversible: bool = True,
) -> list[cobra.Reaction]:
    """Add one porter as facilitated transport of each of its substrates.

    Args:
        model: modified in place; pass a copy or use ``with model:``.
        porter: porter name, used as the heterologous gene identifier. Yeast9
            gene IDs are systematic ORF names; a heterologous gene has no ORF
            name, so the porter name itself is used and is easy to spot.
        substrates: published substrate names (e.g. "D-sorbitol"), translated
            to Yeast9's names via analysis.YEAST9_NAMES.
        reversible: facilitated diffusion carries flux either way; only the
            uptake direction affects growth here.

    Returns:
        The reactions added. Substrates the model cannot represent are skipped
        and reported by ``unrepresentable_substrates``.
    """
    added = []
    for published_name in substrates:
        model_name = analysis.YEAST9_NAMES.get(published_name, published_name)
        if not model_name:
            continue  # cannot be represented; see unrepresentable_substrates

        forms = {
            metabolite.compartment: metabolite
            for metabolite in model.metabolites
            if metabolite.name.strip().lower() == str(model_name).strip().lower()
        }
        outside, inside = forms.get("e"), forms.get("c")
        if outside is None or inside is None:
            continue

        slug = published_name.replace(" ", "_").replace("-", "_")
        reaction = cobra.Reaction(
            f"t_{porter}_{slug}",
            name=f"{published_name} transport via {porter} (facilitated)",
        )
        reaction.add_metabolites({outside: -1, inside: 1})
        reaction.bounds = (-1000.0, 1000.0) if reversible else (0.0, 1000.0)
        # Assigning the rule creates the gene in the model.
        reaction.gene_reaction_rule = porter
        model.add_reactions([reaction])
        added.append(reaction)
    return added


def add_mannitol_porter(model: cobra.Model, porter: str) -> cobra.Reaction:
    """Add mannitol transport, for use with ``add_extracellular_mannitol``.

    Kept separate from ``add_porter`` because mannitol is deliberately mapped
    to "" in ``analysis.YEAST9_NAMES``, so that it can never be modelled by
    accident. Modelling it requires asking for it explicitly, here.
    """
    outside = model.metabolites.get_by_id("s_mannitol_e")
    inside = model.metabolites.get_by_id(CYTOSOLIC_MANNITOL)

    transport = cobra.Reaction(
        f"t_{porter}_D_mannitol",
        name=f"D-mannitol transport via {porter} (facilitated, model extension)",
    )
    transport.add_metabolites({outside: -1, inside: 1})
    transport.bounds = (-1000.0, 1000.0)
    transport.gene_reaction_rule = porter
    model.add_reactions([transport])
    return transport


def unrepresentable_substrates(substrates: list[str]) -> list[str]:
    """Published substrates Yeast9 has no extracellular form for."""
    return [s for s in substrates if not analysis.YEAST9_NAMES.get(s, s)]
