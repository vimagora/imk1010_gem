"""Re-derive every number in README.md from the result CSVs and the data files.

Nothing here trusts the prose: each check recomputes a value from the committed
results and asserts that the README states it. It also makes a few negative
checks -- that no unverified authorship claim has crept back in, that lactose is
not described as modelled when it is not, and that the repaired strain is still
labelled as a diagnostic rather than a real organism.

This exists because a README full of hand-copied numbers goes stale silently.
Run it after any change to the analyses:

    python scripts/run_part_a.py
    python scripts/run_part_b.py
    python scripts/run_part_c.py
    python scripts/verify_readme.py

Exits non-zero if any check fails.
"""

import csv
import hashlib
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

README = (ROOT / "README.md").read_text(encoding="utf-8")
RESULTS = ROOT / "results"
DATA = ROOT / "data"

checks: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    checks.append((name, bool(ok), detail))


def rows(path: pathlib.Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def contains(text: str) -> bool:
    return text in README


# --- model provenance -------------------------------------------------------
prov = json.loads((DATA / "model" / "MODEL_PROVENANCE.json").read_text())
check("commit sha in README", contains(prov["commit_sha"]), prov["commit_sha"][:12])
check("release tag in README", contains(prov["release_tag"]), prov["release_tag"])
check("sha256 in README", contains(prov["sha256"]), prov["sha256"][:12])

model_xml = DATA / "model" / "yeast-GEM.xml"
if model_xml.exists():
    digest = hashlib.sha256(model_xml.read_bytes()).hexdigest()
    check("model file matches recorded sha256", digest == prov["sha256"])
else:
    check("model file present", False, "run scripts/download_model.py")

# --- Part A -----------------------------------------------------------------
panel = rows(RESULTS / "carbon_panel.csv")

for paper, label in [("deValk2022", "de Valk"), ("Nummela2026", "Nummela")]:
    scoreable = [
        r for r in panel
        if r[f"[{paper}] published"] in {"yes", "no"}
        and r[f"[{paper}] IMK1010"] != "not representable"
    ]
    n = len(scoreable)
    always_no = sum(1 for r in scoreable if r[f"[{paper}] published"] == "no")
    for strain, column in [("faithful", "IMK1010"), ("repaired", "IMK1010 repaired")]:
        agree = sum(1 for r in scoreable if r[f"[{paper}] {column}"] == "agrees")
        check(f"Part A {label} {strain} {agree}/{n}", contains(f"{agree}/{n}"),
              f"{agree}/{n}")
    check(f"Part A {label} trivial baseline {always_no}/{n}",
          contains(f"{always_no}/{n}"), f"always-no = {always_no}/{n}")

# The faithful model must be worse than the trivial predictor on both panels,
# which is the README's central claim about the ratios.
for paper in ["deValk2022", "Nummela2026"]:
    scoreable = [
        r for r in panel
        if r[f"[{paper}] published"] in {"yes", "no"}
        and r[f"[{paper}] IMK1010"] != "not representable"
    ]
    faithful = sum(1 for r in scoreable if r[f"[{paper}] IMK1010"] == "agrees")
    always_no = sum(1 for r in scoreable if r[f"[{paper}] published"] == "no")
    check(f"{paper}: faithful worse than trivial predictor", faithful < always_no,
          f"{faithful} < {always_no}")
    repaired = sum(1 for r in scoreable if r[f"[{paper}] IMK1010 repaired"] == "agrees")
    check(f"{paper}: repaired beats trivial predictor", repaired > always_no,
          f"{repaired} > {always_no}")

# Each published panel must contain at least one positive, or the baseline
# argument in the README is misstated.
for paper in ["deValk2022", "Nummela2026"]:
    positives = [r for r in panel if r[f"[{paper}] published"] == "yes"]
    check(f"{paper}: panel contains a published positive", len(positives) >= 1,
          ", ".join(r["substrate"] for r in positives))

# The faithful model's two correct Nummela calls are ethanol and D-galactose.
correct = sorted(
    r["substrate"] for r in panel
    if r["[Nummela2026] IMK1010"] == "agrees"
)
check("Part A: faithful gets ethanol and D-galactose right",
      correct == ["D-galactose", "ethanol"], ", ".join(correct))

# --- survivor genes ---------------------------------------------------------
survivors = rows(RESULTS / "uptake_survivors.csv")
# Keyed by the Yeast9 metabolite name, which is what the CSV records:
# D-sorbitol is stored as D-glucitol and D-xylitol as xylitol.
expected_survivors = {
    "D-glucose": {"CIN10", "HXT17"},
    "D-fructose": {"HXT17"},
    "D-mannose": {"HXT17"},
    "D-glucitol": {"HXT17"},
    "xylitol": {"FPS1"},
    "glycerol": {"GUP1", "GUP2"},
}
for substrate, genes in expected_survivors.items():
    found = {
        g.strip()
        for r in survivors if r["carbon_source"] == substrate
        for g in r["surviving_genes"].split(",") if g.strip()
    }
    check(f"survivors {substrate}", found == genes, f"{sorted(found)}")
    for gene in genes:
        check(f"README names {gene}", contains(gene))

blocked = sorted({
    r["substrate"] for r in survivors
    if r["uptake_blocked_in_IMK1010"] == "yes"
})
check("genuinely blocked substrates", blocked == ["D-galactose", "maltose", "sucrose"],
      ", ".join(blocked))

# --- Part B -----------------------------------------------------------------
glp = rows(RESULTS / "glp_predictions.csv")
scoreable = [r for r in glp if r["verdict"] in {"agrees", "DISAGREES"}]
positives = [r for r in scoreable if r["published"] == "yes"]
negatives = [r for r in scoreable if r["published"] == "no"]
agree = sum(1 for r in scoreable if r["verdict"] == "agrees")
excluded = [r for r in glp if r["verdict"] == "not representable"]

check(f"Part B {agree}/{len(scoreable)}", contains(f"{agree}/{len(scoreable)}"))
check(f"Part B positives {len(positives)}/{len(positives)}",
      contains(f"{len(positives)}/{len(positives)}"))
check(f"Part B negatives {len(negatives)}/{len(negatives)}",
      contains(f"{len(negatives)}/{len(negatives)}"))
check("Part B total pairs stated", contains(str(len(glp))), f"{len(glp)} pairs")
check("Part B all pairs agree", agree == len(scoreable))
check("Part B has published negatives", len(negatives) > 0, f"{len(negatives)}")
check("Part B excluded pairs are mannitol only",
      sorted({r["substrate"] for r in excluded}) == ["D-mannitol"])

# No "not reported" verdicts should survive: they were resolved to negatives.
check("Part B has no unresolved pairs",
      not any(r["verdict"] == "unknown" for r in glp))

# The discriminating asymmetry.
for porter, sugar, expect in [
    ("GlpA", "D-glucose", "no"), ("GlpB", "D-glucose", "no"),
    ("GlpC", "D-glucose", "yes"), ("GlpA", "D-mannose", "no"),
    ("GlpC", "D-mannose", "yes"), ("GlpA", "glycerol", "yes"),
]:
    record = next(r for r in glp if r["porter"] == porter and r["substrate"] == sugar)
    check(f"{porter}/{sugar} published={expect}", record["published"] == expect,
          record["published"])
    check(f"{porter}/{sugar} predicted matches", record["verdict"] == "agrees")

for porter in ["GlpA", "GlpB", "GlpC"]:
    record = next(r for r in glp if r["porter"] == porter
                  and r["substrate"] == "D-galactose")
    check(f"{porter}/D-galactose is a published negative",
          record["published"] == "no")
    record = next(r for r in glp if r["porter"] == porter
                  and r["substrate"] == "ethanol")
    check(f"{porter}/ethanol is the positive control", record["published"] == "yes")

# --- Part C -----------------------------------------------------------------
cello = rows(RESULTS / "cellobiose_panel.csv")
matches = sum(1 for r in cello if r["matches_expectation"].lower() == "yes")
check(f"Part C {matches}/{len(cello)}", contains(f"{matches}/{len(cello)}"))

by_construct = {r["construct"]: r for r in cello}
symport = next(r for r in cello if "symport" in r["construct"])
facilitated = next(r for r in cello if "facilitated" in r["construct"])
penalty = (1 - float(symport["growth_rate"]) / float(facilitated["growth_rate"])) * 100
check(f"Part C symport penalty {penalty:.1f}%", contains(f"{penalty:.1f}%"),
      f"{penalty:.1f}%")
check("Part C symport rate in README",
      contains(f"{float(symport['growth_rate']):.4f}"))
check("Part C facilitated rate in README",
      contains(f"{float(facilitated['growth_rate']):.4f}"))

platform = next(r for r in cello if "no transporter" in r["construct"]
                and "intracellular" in r["construct"])
check("Part C platform strain does not grow", platform["grows"] == "no")
secreted_wt = next(r for r in cello if r["strain"] == "wild type"
                   and "secreted" in r["construct"])
secreted_ko = next(r for r in cello if "IMK1010" in r["strain"]
                   and "secreted" in r["construct"])
check("Part C secreted BGL rescues wild type", secreted_wt["grows"] == "yes")
check("Part C secreted BGL does not rescue IMK1010", secreted_ko["grows"] == "no")

# --- constructs and accessions ----------------------------------------------
constructs = rows(DATA / "cellobiose_constructs.csv")
for accession in ["XP_011393758.1", "MH648002", "NCU00130", "An12g09270"]:
    check(f"accession {accession} in README", contains(accession))
for row in constructs:
    if row["source_organism"].startswith("("):
        continue
    check(f"{row['construct']} cites where it was characterised",
          bool(row["characterised_in"].strip()) and bool(row["doi_or_pmid"].strip()))

# --- gene list ---------------------------------------------------------------
genes = rows(DATA / "imk1010_deleted_genes.csv")
transporters = sum(1 for g in genes if g["category"] == "transporter")
hydrolases = sum(1 for g in genes if g["category"] == "hydrolase")
check("23 transporters", transporters == 23, str(transporters))
check("9 hydrolases", hydrolases == 9, str(hydrolases))
check("32 genes stated in README", contains("32"))
check("29 of 32 map", contains("29 of the 32") or contains("29 of 32"))

# --- negative checks ---------------------------------------------------------
check("no Nummela co-authorship claim",
      not re.search(r"co-author of Nummela", README, re.I))
check("Xu co-authorship retained", bool(re.search(r"co-author of Xu", README, re.I)))
check("lactose not claimed as modelled",
      "not modelled" in README.lower() or "lactose" not in README.lower())
check("culture-format limitation present",
      "culture format" in README.lower())
check("anaerobic attempt reported",
      "anaerobic" in README.lower())
check("solid readout named as the comparison",
      bool(re.search(r"\*\*solid\*\*|solid readout", README)))
check("repaired strain labelled as not real",
      "does not exist in any laboratory" in README
      or "exists in no laboratory" in README)

# --- figures -----------------------------------------------------------------
for stem in ["figure1_part_a_carbon_panel", "figure2_part_b_glp_porters",
             "figure3_part_c_cellobiose"]:
    check(f"{stem}.png exists", (RESULTS / "figures" / f"{stem}.png").exists())
    check(f"{stem} referenced in README", contains(f"{stem}.png"))

# --- report ------------------------------------------------------------------
passed = sum(1 for _, ok, _ in checks if ok)
for name, ok, detail in checks:
    if not ok:
        print(f"FAIL  {name}  {detail}")
print(f"\n{passed}/{len(checks)} checks pass")
sys.exit(0 if passed == len(checks) else 1)
