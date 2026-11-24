"""Part A: can Yeast9 reproduce the published phenotypes of IMK1010?

Three strains are built and each is scored against *both* published sources
separately:

1. IMK1010 -- the 32 published deletions, faithful to the genotype. This is
   the primary result.
2. IMK1010 (repaired) -- additionally deleting the genes Yeast9 lists in
   sugar uptake rules that were never deleted in the real strain. No such
   strain exists in any lab; it isolates the effect of the model's annotation.
3. IMK1010 (repaired, Fps1 conducting inwards) -- as 2, but allowing the
   glycerol channel to work in both directions, as a channel does.

Why score the sources separately: de Valk et al. (2022) saw IMK1010 grow on
2% glycerol in 8 days, while Nummela et al. (2026) saw no growth of the
empty-vector control on 1% glycerol in 6 days. Variant 3 reproduces the first
observation and variant 2 the second. FBA has no concentrations, kinetics or
time, so it cannot reconcile them, and we do not pretend otherwise.

Writes results/carbon_panel.csv and results/uptake_survivors.csv.

Usage:
    python scripts/run_part_a.py
"""

import pandas as pd

from imk1010_gem import analysis, strains
from imk1010_gem.model_io import gene_name_map, load_model

PHENOTYPES = analysis.RESULTS_DIR.parent / "data" / "published_phenotypes.csv"

# (reporting name, strain name in the phenotype file, paper)
SOURCES = [
    ("de Valk 2022 (IMK1010)", "IMK1010", "deValk2022"),
    ("Nummela 2026 (empty vector)", "IMK1010 + empty vector", "Nummela2026"),
]


def main() -> None:
    model = load_model()
    orf_to_name = gene_name_map()

    print(f"Yeast9 loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # --- strain construction -------------------------------------------------
    imk1010, mapping = strains.build_imk1010(model)
    print("\n" + mapping.summary())

    # Repair over every substrate either paper tested, so the repaired strain
    # is not quietly specific to one panel.
    phenotypes = pd.read_csv(PHENOTYPES)
    tested = sorted(set(phenotypes["substrate"]))
    substrates = {name: analysis.YEAST9_NAMES.get(name, name) for name in tested}

    repaired, _, extra = strains.build_repaired_imk1010(
        model, [n for n in substrates.values() if n]
    )

    print("\nGenes Yeast9 leaves holding uptake open in IMK1010:")
    for substrate, orfs in sorted(extra.items()):
        print(f"  {substrate:12s} {', '.join(orf_to_name.get(o, o) for o in orfs)}")

    bidirectional = repaired.copy()
    bidirectional.id = "IMK1010_repaired_Fps1_inward"
    strains.allow_glycerol_channel_uptake(bidirectional)

    variants = {
        "IMK1010": imk1010,
        "IMK1010 repaired": repaired,
        "IMK1010 repaired + Fps1 inward": bidirectional,
    }

    # --- the carbon panel ----------------------------------------------------
    panel = analysis.carbon_panel(model, variants, substrates=substrates)

    for label, published_strain, paper in SOURCES:
        for variant in variants:
            scored = analysis.compare_with_published(
                panel, phenotypes, strain=variant,
                published_strain=published_strain, source_paper=paper,
            )
            panel[f"[{paper}] published"] = scored["published"]
            panel[f"[{paper}] {variant}"] = scored["agreement"]

    print("\n=== Predicted growth ===")
    print(panel[["substrate", "exchange"]
                + [f"{v} grows" for v in variants]].to_string(index=False))

    for label, _, paper in SOURCES:
        cols = [f"[{paper}] {v}" for v in variants]
        sub = panel[panel[f"[{paper}] published"].isin(["yes", "no"])]
        # Substrates Yeast9 cannot represent are excluded from the score and
        # reported separately, so the denominator is not quietly inflated.
        scoreable = sub[sub[f"[{paper}] {list(variants)[0]}"] != "not representable"]
        print(f"\n=== vs {label}: {len(sub)} substrates with a published phenotype ===")
        print(sub[["substrate", f"[{paper}] published"]
                  + [f"{v} grows" for v in variants]].to_string(index=False))
        n_absent = len(sub) - len(scoreable)
        for variant, col in zip(variants, cols):
            print(f"  {variant:32s} {(scoreable[col] == 'agrees').sum()}/{len(scoreable)}")
        if n_absent:
            print(f"  ({n_absent} substrate(s) excluded: not representable in Yeast9)")

    survivors = analysis.uptake_survivor_table(
        model, mapping.orfs, substrates=substrates, orf_to_name=orf_to_name
    )
    for frame, name in [(panel, "carbon_panel.csv"), (survivors, "uptake_survivors.csv")]:
        print(f"\nwrote {analysis.write_csv(frame, name)}")


if __name__ == "__main__":
    main()
