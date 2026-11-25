"""Part B: do the A. niger Glp porters restore growth of IMK1010 in silico?

Background strain: IMK1010 *repaired*. In the genotype-faithful strain Yeast9
already predicts growth on all of these substrates (Part A), so no porter could
show an effect. The repaired strain is the one that reproduces Nummela's
empty-vector control, which grew only on ethanol.

Both published positives and published negatives are scored. The negatives are
the informative half: GlpA and GlpB were tested on the hexoses and did not
restore growth, while GlpC did. A model that simply said "yes" to every added
transporter would get all the positives right and every negative wrong.

Comparison is against the SOLID-medium readout; see the README limitations for
why, and for what the liquid-culture results show.

Writes results/glp_predictions.csv.

Usage:
    python scripts/run_part_b.py
"""

import pandas as pd

from imk1010_gem import analysis, porters, strains


def main() -> None:
    from imk1010_gem.model_io import load_model

    model = load_model()
    published = pd.read_csv(porters.GLP_CSV)
    substrates = sorted(set(published["substrate"]))

    representable = [
        analysis.YEAST9_NAMES.get(s, s)
        for s in substrates
        if analysis.YEAST9_NAMES.get(s, s)
    ]
    repaired, _, _ = strains.build_repaired_imk1010(model, representable)
    print(f"background strain: {repaired.id}")

    missing = porters.unrepresentable_substrates(substrates)
    if missing:
        print(f"substrates Yeast9 cannot represent: {', '.join(missing)}")

    # Substrates each porter was published to carry: what we add to the model.
    specificity = porters.read_glp_specificity()

    records = []
    for porter in sorted(specificity):
        carried = specificity[porter]
        print(f"\n{porter}: adding transport of {', '.join(sorted(carried))}")

        for substrate in substrates:
            model_name = analysis.YEAST9_NAMES.get(substrate, substrate)
            observed = published[
                (published["porter"] == porter) & (published["substrate"] == substrate)
            ]["transported"].iloc[0]

            if not model_name:
                predicted, rate, verdict = "not representable", 0.0, "not representable"
            else:
                with repaired as strain:
                    porters.add_porter(strain, porter, carried)
                    rate = analysis.growth_on(strain, model_name)
                predicted = "yes" if analysis.grows(rate) else "no"
                if observed in {"yes", "no"}:
                    verdict = "agrees" if predicted == observed else "DISAGREES"
                else:
                    verdict = "unknown"

            records.append({
                "porter": porter, "substrate": substrate,
                "published": observed, "predicted": predicted,
                "growth_rate": round(rate, 6), "verdict": verdict,
            })
            print(f"   {substrate:12s} published={observed:4s} "
                  f"predicted={predicted:17s} {verdict}")

    table = pd.DataFrame.from_records(records)
    scoreable = table[table["verdict"].isin(["agrees", "DISAGREES"])]
    n_agree = (scoreable["verdict"] == "agrees").sum()

    positives = scoreable[scoreable["published"] == "yes"]
    negatives = scoreable[scoreable["published"] == "no"]

    print(f"\n=== Agreement with Nummela et al. (2026), solid medium ===")
    print(f"  overall                 {n_agree}/{len(scoreable)} porter-substrate pairs")
    print(f"    published positives   {(positives['verdict'] == 'agrees').sum()}/{len(positives)}")
    print(f"    published negatives   {(negatives['verdict'] == 'agrees').sum()}/{len(negatives)}")

    wrong = scoreable[scoreable["verdict"] == "DISAGREES"]
    if len(wrong):
        print("  disagrees on: "
              + ", ".join(f"{r.porter}/{r.substrate}" for r in wrong.itertuples()))

    excluded = table[table["verdict"] == "not representable"]
    if len(excluded):
        print(f"  {len(excluded)} pairs excluded, not representable in Yeast9: "
              f"{', '.join(sorted(set(excluded['substrate'])))}")

    print(f"\nwrote {analysis.write_csv(table, 'glp_predictions.csv')}")


if __name__ == "__main__":
    main()
