"""Part A: can Yeast9 reproduce the published phenotypes of IMK1010?

Builds two strains and compares both with de Valk et al. (2022):

1. IMK1010 -- the 32 published deletions, faithful to the genotype.
2. IMK1010 (repaired) -- additionally deleting the genes Yeast9 lists in
   sugar uptake rules that were never deleted in the real strain. No such
   strain exists; it isolates the effect of the model's gene annotation.

Writes results/carbon_panel.csv and results/uptake_survivors.csv.

Usage:
    python scripts/run_part_a.py
"""

import pandas as pd

from imk1010_gem import analysis, strains
from imk1010_gem.model_io import gene_name_map, load_model

PHENOTYPES = analysis.RESULTS_DIR.parent / "data" / "devalk2022_phenotypes.csv"


def main() -> None:
    model = load_model()
    orf_to_name = gene_name_map()

    print(f"Yeast9 loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # --- strain construction -------------------------------------------------
    imk1010, mapping = strains.build_imk1010(model)
    print("\n" + mapping.summary())

    substrates = list(analysis.DEVALK_SUBSTRATES.values())
    repaired, _, extra = strains.build_repaired_imk1010(model, substrates)

    print("\nGenes Yeast9 leaves holding sugar uptake open in IMK1010:")
    for substrate, orfs in extra.items():
        names = ", ".join(orf_to_name.get(o, o) for o in orfs)
        print(f"  {substrate:12s} {names}")
    if not extra:
        print("  (none)")

    # A third variant: the repaired strain, but with the glycerol channel
    # allowed to conduct inwards as a real channel would. See
    # strains.allow_glycerol_channel_uptake for the reasoning.
    bidirectional = repaired.copy()
    bidirectional.id = "IMK1010_repaired_bidirectional_Fps1"
    channel = strains.allow_glycerol_channel_uptake(bidirectional)
    print(f"\nglycerol channel {channel.id} bounds in third variant: {channel.bounds}")

    # --- the carbon panel ----------------------------------------------------
    panel = analysis.carbon_panel(
        model,
        {
            "IMK1010": imk1010,
            "IMK1010 repaired": repaired,
            "IMK1010 repaired + Fps1 uptake": bidirectional,
        },
    )
    phenotypes = pd.read_csv(PHENOTYPES)
    panel = analysis.compare_with_published(panel, phenotypes, strain="IMK1010")

    # Score each variant against the same published phenotypes.
    variants = ["IMK1010", "IMK1010 repaired", "IMK1010 repaired + Fps1 uptake"]
    for variant in variants[1:]:
        scored = analysis.compare_with_published(
            panel.drop(columns=["published", "agreement"], errors="ignore"),
            phenotypes,
            strain=variant,
        )
        panel[f"{variant} agreement"] = scored["agreement"]
    panel = panel.rename(columns={"agreement": "IMK1010 agreement"})

    columns = (
        ["substrate", "exchange", "published", "wild type grows"]
        + [c for v in variants for c in (f"{v} grows", f"{v} agreement")]
        + ["wild type growth"] + [f"{v} growth" for v in variants]
    )
    panel = panel[columns]

    print("\n=== Carbon panel: predicted vs published ===")
    print(panel[["substrate", "published"]
                + [c for v in variants for c in (f"{v} grows", f"{v} agreement")]
                ].to_string(index=False))

    survivors = analysis.uptake_survivor_table(
        model, mapping.orfs, orf_to_name=orf_to_name
    )

    for frame, name in [(panel, "carbon_panel.csv"), (survivors, "uptake_survivors.csv")]:
        print(f"\nwrote {analysis.write_csv(frame, name)}")

    known = panel[panel["IMK1010 agreement"] != "unknown"]
    if len(known):
        print(f"\nAgreement with de Valk et al. (2022), "
              f"{len(known)} substrates with a published phenotype:")
        for variant in variants:
            n_agree = (known[f"{variant} agreement"] == "agrees").sum()
            print(f"  {variant:34s} {n_agree}/{len(known)}")


if __name__ == "__main__":
    main()
