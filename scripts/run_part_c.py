"""Part C: is transport the limiting step for cellobiose utilisation?

S. cerevisiae lacks both a cellobiose transporter and a beta-glucosidase, so
two things are missing. The platform strain IMK1010 + NcBGL supplies only the
enzyme, which isolates the transporter as the variable.

Strains compared:
  wild type                      - context
  IMK1010 (repaired)             - the background; grows on none of these
  IMK1010 + NcBGL                - enzyme only, no transporter
  IMK1010 + NcBGL + CtA          - enzyme and transporter

Plus the localisation control: a secreted beta-glucosidase in each background.

Writes results/cellobiose_panel.csv.

Usage:
    python scripts/run_part_c.py
"""

import pandas as pd

from imk1010_gem import analysis, cellobiose, strains
from imk1010_gem.model_io import load_model

SUBSTRATES = ["D-glucose", "D-fructose", "D-mannose", "D-glucitol", "xylitol", "glycerol"]


def main() -> None:
    model = load_model()
    constructs = cellobiose.read_constructs()
    print("Heterologous constructs:")
    for row in constructs:
        print(f"  {row['construct']:12s} {row['type']:32s} "
              f"{row['source_organism']} ({row['gene_or_locus']})")

    repaired, _, _ = strains.build_repaired_imk1010(model, SUBSTRATES)
    records = []

    def record(strain_name, construct, rate, expectation):
        grew = analysis.grows(rate)
        records.append({
            "strain": strain_name, "construct": construct,
            "growth_rate": round(rate, 6), "grows": "yes" if grew else "no",
            "expected": expectation,
            "matches_expectation": "yes" if (grew == (expectation == "yes")) else "NO",
        })
        print(f"  {strain_name:22s} {construct:34s} {rate:.4f}  "
              f"{'GROWS' if grew else 'no growth'}")

    print("\n=== Growth on cellobiose ===")

    # 1. Backgrounds without any construct: cellobiose does not exist for them.
    for name, base in [("wild type", model), ("IMK1010 repaired", repaired)]:
        with base as strain:
            cellobiose.add_cellobiose_pools(strain)
            rate = analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME)
        record(name, "cellobiose available, no enzyme", rate, "no")

    # 2. The platform strain: enzyme but no transporter. Must not grow.
    with repaired as strain:
        _, inside = cellobiose.add_cellobiose_pools(strain)
        cellobiose.add_intracellular_bgl(strain, inside=inside)
        rate = analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME)
    record("IMK1010 + NcBGL", "intracellular BGL, no transporter", rate, "no")

    # 3. Enzyme plus transporter, under both mechanistic assumptions.
    for symport in [True, False]:
        label = "H+ symport" if symport else "facilitated"
        with repaired as strain:
            _, inside = cellobiose.add_cellobiose_pools(strain)
            cellobiose.add_intracellular_bgl(strain, inside=inside)
            cellobiose.add_cellobiose_transporter(strain, "CtA", proton_symport=symport)
            rate = analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME)
        record("IMK1010 + NcBGL + CtA", f"intracellular BGL + CtA ({label})", rate, "yes")

    # 4. Localisation control: a secreted BGL instead of an intracellular one.
    print("\n=== Localisation control: secreted beta-glucosidase ===")
    for name, base, expectation in [("wild type", model, "yes"),
                                     ("IMK1010 repaired", repaired, "no")]:
        with base as strain:
            cellobiose.add_cellobiose_pools(strain)
            cellobiose.add_secreted_bgl(strain)
            rate = analysis.growth_on(strain, cellobiose.CELLOBIOSE_NAME)
        record(name, "secreted BGL, no transporter", rate, expectation)

    table = pd.DataFrame.from_records(records)
    failures = (table["matches_expectation"] == "NO").sum()
    print(f"\n{len(table) - failures}/{len(table)} results match expectation")

    symport = table[table["construct"].str.contains("H. symport", regex=True)]
    facilitated = table[table["construct"].str.contains("facilitated")]
    if len(symport) and len(facilitated):
        s, f = symport["growth_rate"].iloc[0], facilitated["growth_rate"].iloc[0]
        print(f"proton symport costs {(1 - s / f) * 100:.1f}% of growth rate "
              f"({s:.4f} vs {f:.4f}); per hexose that is half the cost of "
              f"proton-coupled glucose uptake, since one event carries two")

    print(f"\nwrote {analysis.write_csv(table, 'cellobiose_panel.csv')}")


if __name__ == "__main__":
    main()
