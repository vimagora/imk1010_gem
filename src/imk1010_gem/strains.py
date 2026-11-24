"""Building IMK1010 in silico.

Two strains are built here, and the distinction matters:

``build_imk1010``
    The *genotype-faithful* strain: exactly the 32 genes de Valk et al. (2022)
    deleted, knocked out in Yeast9. This is the honest reconstruction and the
    one we report as the primary result.

``build_repaired_imk1010``
    The same strain plus any *additional* genes that Yeast9 lists in a sugar
    transport gene rule but that were never deleted in the real strain. This
    strain does not exist in any lab. It is built only to show what the model
    would predict if its transport gene rules matched the strain's phenotype.

The second model exists because of a real finding. Yeast9's gene rules for
sugar uptake are long OR clauses, and an OR clause survives if a single gene
survives. Yeast9's glucose transport reaction (r_1166) lists twenty genes
including HXT17 and CIN10, neither of which de Valk deleted -- so deleting all
32 published genes leaves glucose uptake completely open, and the model
predicts growth where the experiment shows none.

Crucially, the extra genes are *derived from the model* by
``surviving_uptake_genes`` rather than chosen by hand, so the repair is
reproducible and cannot be quietly tuned to make results look better.
"""

from __future__ import annotations

import csv
import pathlib
from dataclasses import dataclass, field

import cobra

from .model_io import common_name_to_orf

DELETED_GENES_CSV = (
    pathlib.Path(__file__).resolve().parents[2] / "data" / "imk1010_deleted_genes.csv"
)


@dataclass
class GeneMapping:
    """The result of mapping published gene names onto model gene IDs.

    ``unmapped`` is reported rather than silently dropped: a gene that is
    missing from the model is a difference between the model's reference
    genome (S288C) and the strain's background (CEN.PK), which is a finding.
    """

    mapped: dict[str, str] = field(default_factory=dict)  # common name -> ORF
    unmapped: list[str] = field(default_factory=list)

    @property
    def orfs(self) -> list[str]:
        return sorted(self.mapped.values())

    def summary(self) -> str:
        lines = [f"{len(self.mapped)} of {len(self.mapped) + len(self.unmapped)} genes map to Yeast9"]
        if self.unmapped:
            lines.append(f"  not in the model: {', '.join(self.unmapped)}")
        return "\n".join(lines)


def read_deleted_genes(path: pathlib.Path | str = DELETED_GENES_CSV) -> list[dict[str, str]]:
    """Read the published deleted-gene list from the data file."""
    with pathlib.Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def map_genes_to_model(model: cobra.Model, gene_names: list[str]) -> GeneMapping:
    """Map common gene names onto Yeast9 gene IDs, reporting failures."""
    name_to_orf = common_name_to_orf()
    model_genes = {gene.id for gene in model.genes}

    mapping = GeneMapping()
    for name in gene_names:
        orf = name_to_orf.get(name.upper())
        if orf is not None and orf in model_genes:
            mapping.mapped[name] = orf
        else:
            mapping.unmapped.append(name)
    return mapping


def transport_reactions(
    model: cobra.Model, metabolite_name: str
) -> list[cobra.Reaction]:
    """Reactions that move a metabolite between the extracellular space and cytosol.

    Found by looking for reactions containing *both* the [e] and [c] forms of
    the metabolite, so we do not depend on reaction naming.
    """
    forms = {
        metabolite.compartment: metabolite
        for metabolite in model.metabolites
        if metabolite.name.strip().lower() == metabolite_name.strip().lower()
    }
    outside, inside = forms.get("e"), forms.get("c")
    if outside is None or inside is None:
        return []
    return sorted(set(outside.reactions) & set(inside.reactions), key=lambda r: r.id)


def can_import(reaction: cobra.Reaction, metabolite_name: str) -> bool:
    """Whether this reaction's bounds currently permit *net uptake*.

    Direction matters. Yeast9's glycerol channel (r_1172) is written
    cytosol -> extracellular and is irreversible, so it can only export; it
    must not be counted as an uptake route.
    """
    inside = outside = None
    for metabolite, coefficient in reaction.metabolites.items():
        if metabolite.name.strip().lower() != metabolite_name.strip().lower():
            continue
        if metabolite.compartment == "c":
            inside = coefficient
        elif metabolite.compartment == "e":
            outside = coefficient
    if inside is None or outside is None:
        return False
    # Written as uptake (consumes [e], produces [c]): needs forward flux.
    if outside < 0 < inside:
        return reaction.upper_bound > 0
    # Written as export: uptake requires running it in reverse.
    if inside < 0 < outside:
        return reaction.lower_bound < 0
    return False


def surviving_uptake_genes(
    model: cobra.Model, metabolite_name: str, deleted_orfs: list[str]
) -> dict[str, list[str]]:
    """Genes that keep a substrate's uptake open despite the deletions.

    Returns ``{reaction_id: [surviving gene ORFs]}`` for each uptake reaction
    that is still active after the deletions. An empty dict means uptake of
    this substrate is genuinely blocked in the knockout strain.

    This is the function that *derives* the repair set, so that the second
    model is a consequence of the model's own annotation rather than our
    choice.
    """
    deleted = set(deleted_orfs)
    survivors: dict[str, list[str]] = {}

    with model:
        for orf in deleted:
            if orf in {gene.id for gene in model.genes}:
                model.genes.get_by_id(orf).knock_out()

        for reaction in transport_reactions(model, metabolite_name):
            if not can_import(reaction, metabolite_name):
                continue
            remaining = sorted(
                gene.id for gene in reaction.genes if gene.id not in deleted
            )
            # A reaction with no gene rule at all stays open too: that is
            # passive diffusion in the model, and worth flagging explicitly.
            survivors[reaction.id] = remaining

    return survivors


# Yeast9's glycerol channel (Fps1/Aqy3). It is written cytosol -> extracellular
# and irreversible, so in the model glycerol can only leave the cell this way.
GLYCEROL_CHANNEL = "r_1172"


def allow_glycerol_channel_uptake(model: cobra.Model) -> cobra.Reaction:
    """Make the glycerol channel bidirectional.

    Fps1 is a channel, and channels conduct in both directions; which way
    glycerol actually moves is set by the concentration gradient, not by the
    protein. Yeast9 nonetheless constrains r_1172 to export only, which means
    the model has *no* passive glycerol uptake route at all.

    KEPT FOR THE RECORD, BUT NOT NEEDED BY THE DATA.

    We first reached for this because de Valk et al. appeared to show IMK1010
    growing on glycerol despite lacking STL1, which would have required a
    passive uptake route the model does not have. On a closer reading, their
    ethanol/glycerol plate contained *both* carbon sources, and ethanol alone
    accounts for the growth -- so there is no published evidence of
    glycerol-only growth, and Nummela et al. (2026) report none for the
    empty-vector control. No observation requires this edit.

    It stays here because it is a real limitation of Yeast9 worth naming (a
    channel modelled as one-way), and because the hypothesis and its
    withdrawal belong in the record rather than being quietly deleted.
    """
    reaction = model.reactions.get_by_id(GLYCEROL_CHANNEL)
    reaction.lower_bound = -abs(reaction.upper_bound)
    return reaction


def _knock_out(model: cobra.Model, orfs: list[str]) -> None:
    present = {gene.id for gene in model.genes}
    for orf in orfs:
        if orf in present:
            model.genes.get_by_id(orf).knock_out()


def build_imk1010(
    model: cobra.Model,
    path: pathlib.Path | str = DELETED_GENES_CSV,
    express: tuple[str, ...] = (),
) -> tuple[cobra.Model, GeneMapping]:
    """Return a copy of the model with the 32 published genes knocked out.

    Works on a copy so the wild-type model stays available for comparison in
    the same session.

    Args:
        express: common names of deleted genes to leave *functional*, for
            modelling a complementation strain. Nummela et al. used
            HXT5-expressing IMK1010 as a positive control, which is this
            deletion background carrying HXT5 on a plasmid; in flux terms
            that is the same as not deleting it.
    """
    names = [row["common_name"] for row in read_deleted_genes(path)]
    mapping = map_genes_to_model(model, names)

    restored = {name.upper() for name in express}
    to_delete = [orf for name, orf in mapping.mapped.items() if name.upper() not in restored]

    strain = model.copy()
    strain.id = "IMK1010" + ("+" + "+".join(sorted(restored)) if restored else "")
    _knock_out(strain, to_delete)
    return strain, mapping


def build_repaired_imk1010(
    model: cobra.Model,
    substrates: list[str],
    path: pathlib.Path | str = DELETED_GENES_CSV,
    express: tuple[str, ...] = (),
) -> tuple[cobra.Model, GeneMapping, dict[str, list[str]]]:
    """IMK1010 plus the genes Yeast9 leaves holding sugar uptake open.

    Args:
        model: wild-type model.
        substrates: substrate names whose uptake routes should be closed.

    Returns:
        ``(strain, mapping, extra)`` where ``extra`` maps each substrate to the
        additional gene ORFs that had to be deleted. ``extra`` is the evidence
        for the finding and belongs in the results table.
    """
    names = [row["common_name"] for row in read_deleted_genes(path)]
    mapping = map_genes_to_model(model, names)

    extra: dict[str, list[str]] = {}
    for substrate in substrates:
        survivors = surviving_uptake_genes(model, substrate, mapping.orfs)
        genes = sorted({orf for orfs in survivors.values() for orf in orfs})
        if genes:
            extra[substrate] = genes

    restored = {name.upper() for name in express}
    name_to_orf = common_name_to_orf()
    restored_orfs = {name_to_orf.get(n) for n in restored}

    everything = set(mapping.orfs) | {o for orfs in extra.values() for o in orfs}
    to_delete = sorted(everything - restored_orfs)

    strain = model.copy()
    strain.id = "IMK1010_repaired" + ("+" + "+".join(sorted(restored)) if restored else "")
    _knock_out(strain, to_delete)
    return strain, mapping, extra
