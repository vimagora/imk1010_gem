"""Part B: do the A. niger Glp porters restore growth of IMK1010 in silico?

Background strain: IMK1010 *repaired*. That is the deliberate choice. In the
genotype-faithful strain Yeast9 already predicts growth on every one of these
substrates (see Part A), so a porter could not possibly show an effect -- there
would be nothing to restore. The repaired strain is the one that reproduces
Nummela's empty-vector control, which is the background the porters were
actually expressed in.

For each porter we add facilitated transport of its published substrates and
ask whether IMK1010 grows. Compared against Nummela et al. (2026) Fig. 4.

Writes results/glp_predictions.csv.

Usage:
    python scripts/run_part_b.py
"""

import pandas as pd

from imk1010_gem import analysis, porters, strains
from imk1010_gem.model_io import load_model

GLP = analysis.RESULTS_DIR.parent / "data" / "nummela2026_glp.csv"


def main() -> None:
    model = load_model()
    specificity = porters.read_glp_specificity()
    published = pd.read_csv(GLP)

    all_substrates = sorted(set(published["substrate"]))
    representable = [
        analysis.YEAST9_NAMES.get(s, s)
        for s in all_substrates
        if analysis.YEAST9_NAMES.get(s, s)
    ]

    repaired, _, _ = strains.build_repaired_imk1010(model, representable)
    print(f"background strain: {repaired.id}")

    missing = porters.unrepresentable_substrates(all_substrates)
    if missing:
        print(f"substrates Yeast9 cannot represent: {', '.join(missing)}")

    records = []
    for porter, substrates in sorted(specificity.items()):
        print(f"\n{porter}: published substrates {', '.join(substrates)}")

        for substrate in all_substrates:
            model_name = analysis.YEAST9_NAMES.get(substrate, substrate)
            published_grows = (
                "yes" if substrate in substrates else "not reported"
            )

            if not model_name:
                records.append({
                    "porter": porter, "substrate": substrate,
                    "published": published_grows,
                    "predicted": "not representable",
                    "verdict": "model cannot represent this substrate",
                })
                continue

            # Add only this porter, then test this substrate.
            with repaired as strain:
                porters.add_porter(strain, porter, substrates)
                rate = analysis.growth_on(strain, model_name)

            predicted = "yes" if analysis.grows(rate) else "no"
            if published_grows == "yes":
                verdict = "agrees" if predicted == "yes" else "DISAGREES"
            else:
                # Not reported is not the same as "does not transport": absence
                # of a published result is not evidence either way.
                verdict = "no published result" if predicted == "no" else \
                          "model predicts growth, not reported in the paper"

            records.append({
                "porter": porter, "substrate": substrate,
                "published": published_grows, "predicted": predicted,
                "growth_rate": round(rate, 6), "verdict": verdict,
            })
            print(f"   {substrate:12s} published={published_grows:12s} "
                  f"predicted={predicted:3s} ({rate:.4f})  {verdict}")

    table = pd.DataFrame.from_records(records)

    claimed = table[table["published"] == "yes"]
    scoreable = claimed[claimed["predicted"] != "not representable"]
    n_agree = (scoreable["verdict"] == "agrees").sum()
    print(f"\n=== Of the {len(claimed)} published porter-substrate pairs ===")
    print(f"  representable in Yeast9 : {len(scoreable)}")
    print(f"  reproduced by the model : {n_agree}/{len(scoreable)}")
    if len(claimed) > len(scoreable):
        print(f"  not representable       : {len(claimed) - len(scoreable)} "
              f"({', '.join(sorted(set(claimed[claimed['predicted'] == 'not representable']['substrate'])))})")

    # --- the mannitol gap ----------------------------------------------------
    # All three porters restored growth on mannitol in the paper, but Yeast9
    # has no extracellular mannitol and no exchange, so it cannot represent
    # that at all. The catabolic step IS present (r_4183). Extending the model
    # with the uptake side separates "missing transporter" from "missing
    # pathway". This is an extension beyond Yeast9, reported as such.
    print("\n=== Mannitol: is the gap the transporter or the pathway? ===")
    extended = repaired.copy()
    extended.id = "IMK1010_repaired_plus_mannitol_uptake"
    porters.add_extracellular_mannitol(extended)

    for porter in sorted(specificity):
        with extended as strain:
            porters.add_mannitol_porter(strain, porter)
            rate = analysis.growth_on(strain, "Mannitol")
        verdict = "grows" if analysis.grows(rate) else "still no growth"
        print(f"   {porter}: {verdict} ({rate:.4f})")

    print("\n  Interpretation: Yeast9 already carries the catabolic step")
    print("  (r_4183, mannitol + NAD <=> D-fructose + NADH, DSF1 or MAN2).")
    print("  If growth appears once uptake is added, the gap in Yeast9 is the")
    print("  transport side only, not the pathway.")

    print(f"\nwrote {analysis.write_csv(table, 'glp_predictions.csv')}")


if __name__ == "__main__":
    main()
