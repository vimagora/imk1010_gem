"""Part A: can Yeast9 reproduce the published phenotypes of IMK1010?

Strains built, and why:

1. IMK1010 -- the 32 published deletions, faithful to the genotype. The
   primary result.
2. IMK1010 (repaired) -- additionally deleting the genes Yeast9 lists in
   sugar uptake rules that were never deleted in the real strain. No such
   strain exists in any lab; it isolates the effect of the model's
   annotation. The extra genes are derived from the model's own gene rules,
   not chosen by hand.
3. IMK1010 + HXT5 -- Nummela et al.'s positive control, modelled by leaving
   HXT5 functional in the repaired background. A genuine positive control:
   if this does not restore hexose growth, our repair has broken something.

Scored against each published source separately, never blended.

Writes results/carbon_panel.csv and results/uptake_survivors.csv.

Usage:
    python scripts/run_part_a.py
"""

import pandas as pd

from imk1010_gem import analysis, strains
from imk1010_gem.model_io import gene_name_map, load_model

PHENOTYPES = analysis.RESULTS_DIR.parent / "data" / "published_phenotypes.csv"

SOURCES = [
    ("de Valk 2022 (IMK1010)", "IMK1010", "deValk2022"),
    ("Nummela 2026 (empty vector)", "IMK1010 + empty vector", "Nummela2026"),
]


def main() -> None:
    model = load_model()
    orf_to_name = gene_name_map()
    phenotypes = pd.read_csv(PHENOTYPES)

    print(f"Yeast9 loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    imk1010, mapping = strains.build_imk1010(model)
    print("\n" + mapping.summary())

    # Every substrate either paper tested, so the repair is not tuned to one panel.
    tested = sorted(set(phenotypes["substrate"]))
    substrates = {name: analysis.YEAST9_NAMES.get(name, name) for name in tested}
    flat = sorted({n for v in substrates.values() if v
                   for n in ([v] if isinstance(v, str) else v)})

    repaired, _, extra = strains.build_repaired_imk1010(model, flat)
    positive, _, _ = strains.build_repaired_imk1010(model, flat, express=("HXT5",))

    print("\nGenes Yeast9 leaves holding uptake open in IMK1010:")
    for substrate, orfs in sorted(extra.items()):
        print(f"  {substrate:12s} {', '.join(orf_to_name.get(o, o) for o in orfs)}")

    variants = {
        "IMK1010": imk1010,
        "IMK1010 repaired": repaired,
        "IMK1010 repaired + HXT5": positive,
    }
    # The HXT5 strain is a different strain from the empty-vector control, so
    # it is reported separately rather than scored against the control's rows.
    scored_variants = ["IMK1010", "IMK1010 repaired"]

    panel = analysis.carbon_panel(model, variants, substrates=substrates)

    for _, published_strain, paper in SOURCES:
        for variant in scored_variants:
            scored = analysis.compare_with_published(
                panel, phenotypes, strain=variant,
                published_strain=published_strain, source_paper=paper,
            )
            panel[f"[{paper}] published"] = scored["published"]
            panel[f"[{paper}] {variant}"] = scored["agreement"]

    print("\n=== Predicted growth ===")
    print(panel[["substrate", "exchange"] + [f"{v} grows" for v in variants]]
          .to_string(index=False))

    for label, _, paper in SOURCES:
        sub = panel[panel[f"[{paper}] published"].isin(["yes", "no"])]
        scoreable = sub[sub[f"[{paper}] {list(variants)[0]}"] != "not representable"]
        print(f"\n=== vs {label} ===")
        print(sub[["substrate", f"[{paper}] published"]
                  + [f"{v} grows" for v in scored_variants]].to_string(index=False))
        for variant in scored_variants:
            n = (scoreable[f"[{paper}] {variant}"] == "agrees").sum()
            print(f"  {variant:28s} {n}/{len(scoreable)}")
        if len(sub) > len(scoreable):
            print(f"  ({len(sub) - len(scoreable)} excluded: not representable in Yeast9)")

    # Positive control: HXT5 is a hexose transporter, so restoring it should
    # restore hexose growth and nothing else. If it restored everything, our
    # repair would be suspect.
    restored = [row["substrate"] for _, row in panel.iterrows()
                if row["IMK1010 repaired + HXT5 grows"] == "yes"
                and row["IMK1010 repaired grows"] == "no"]
    print(f"\n=== Positive control: expressing HXT5 in the repaired strain ===")
    print(f"  restores growth on: {', '.join(restored) or '(nothing)'}")

    survivors = analysis.uptake_survivor_table(
        model, mapping.orfs, substrates=substrates, orf_to_name=orf_to_name
    )
    for frame, name in [(panel, "carbon_panel.csv"), (survivors, "uptake_survivors.csv")]:
        print(f"\nwrote {analysis.write_csv(frame, name)}")


if __name__ == "__main__":
    main()
