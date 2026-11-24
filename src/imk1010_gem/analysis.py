"""Running the comparisons and writing the results tables.

Everything here is deliberately thin: set a medium, solve, record the number.
The judgement lives in `media.py` (what the medium is) and `strains.py` (what
the strain is), so that this module stays easy to read against the results.
"""

from __future__ import annotations

import pathlib
import warnings

import cobra
import pandas as pd

from . import media, strains

RESULTS_DIR = pathlib.Path(__file__).resolve().parents[2] / "results"

# A flux below this is numerical noise from the LP solver, not growth.
GROWTH_THRESHOLD = 1e-6

# The carbon sources de Valk et al. (2022) tested, with the names Yeast9 uses.
DEVALK_SUBSTRATES = {
    "D-glucose": "D-glucose",
    "D-fructose": "D-fructose",
    "D-galactose": "D-galactose",
    "sucrose": "sucrose",
    "maltose": "maltose",
    "ethanol": "ethanol",
    "glycerol": "glycerol",
}

# Published substrate name -> the name Yeast9 stores it under. Yeast9 uses some
# synonyms, and one substrate is simply absent, so this mapping is explicit
# rather than assumed. An empty string means the model cannot represent it.
YEAST9_NAMES = {
    "D-glucose": "D-glucose",
    "D-fructose": "D-fructose",
    "D-galactose": "D-galactose",
    "D-mannose": "D-mannose",
    "sucrose": "sucrose",
    "maltose": "maltose",
    "ethanol": "ethanol",
    "glycerol": "glycerol",
    "D-sorbitol": "D-glucitol",  # same molecule, different accepted name
    "D-xylitol": "xylitol",
    # Mannitol exists only as a cytosolic metabolite (s_3805) with no
    # extracellular form and no exchange reaction, so Yeast9 cannot take it up
    # at all. Treated as unrepresentable and reported as such.
    "D-mannitol": "",
    # de Valk et al. tested a single plate containing BOTH ethanol and
    # glycerol, so growth on it cannot be attributed to either alone. Modelled
    # as a medium with both carbon sources open.
    "ethanol + glycerol": ["ethanol", "glycerol"],
}


def growth_on(
    model: cobra.Model,
    substrate: str | list[str],
    uptake: float = media.DEFAULT_UPTAKE,
    per_carbon: bool = False,
) -> float:
    """Predicted growth rate (h^-1) on ``substrate`` as the only carbon source.

    ``substrate`` may be a list of names, for a medium offering more than one
    carbon source (de Valk's combined ethanol + glycerol plate).

    Returns 0.0 for an infeasible problem, which is what "cannot grow" looks
    like to the solver.
    """
    with model:
        names = [substrate] if isinstance(substrate, str) else list(substrate)
        media.set_carbon_sources(model, names, uptake=uptake, per_carbon=per_carbon)
        # An infeasible problem is an expected outcome here -- it is what "this
        # strain cannot grow on this carbon source" looks like to the solver.
        # We check the status ourselves, so cobrapy's warning is just noise.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Solver status is 'infeasible'")
            solution = model.optimize()
        if solution.status != "optimal" or solution.objective_value is None:
            return 0.0
        return max(0.0, solution.objective_value)


def grows(growth_rate: float) -> bool:
    return growth_rate > GROWTH_THRESHOLD


def carbon_panel(
    wild_type: cobra.Model,
    models: dict[str, cobra.Model],
    substrates: dict[str, str] | None = None,
    uptake: float = media.DEFAULT_UPTAKE,
) -> pd.DataFrame:
    """Predicted growth of each strain on each carbon source.

    Args:
        wild_type: the unmodified model, reported as the reference.
        models: other strains to test, keyed by the name to report them under.
        substrates: mapping of published substrate name -> Yeast9 metabolite
            name, defaulting to de Valk's panel.
        uptake: uptake bound applied to whichever carbon source is open.
    """
    substrates = substrates or DEVALK_SUBSTRATES
    all_models = {"wild type": wild_type, **models}

    records = []
    for published_name, model_name in substrates.items():
        if not model_name:
            # Yeast9 cannot represent this substrate at all (e.g. mannitol has
            # no extracellular form). Record it rather than dropping it.
            records.append({
                "substrate": published_name,
                "model_metabolite": "(absent from Yeast9)",
                "exchange": "",
                **{f"{n} growth": 0.0 for n in all_models},
                **{f"{n} grows": "not representable" for n in all_models},
            })
            continue
        row: dict[str, object] = {
            "substrate": published_name,
            "model_metabolite": model_name,
            "exchange": ", ".join(
                media.find_exchange(wild_type, n).id
                for n in ([model_name] if isinstance(model_name, str) else model_name)
            ),
        }
        for strain_name, model in all_models.items():
            rate = growth_on(model, model_name, uptake=uptake)
            row[f"{strain_name} growth"] = round(rate, 6)
            row[f"{strain_name} grows"] = "yes" if grows(rate) else "no"
        records.append(row)
    return pd.DataFrame.from_records(records)


def compare_with_published(
    panel: pd.DataFrame,
    phenotypes: pd.DataFrame,
    strain: str = "IMK1010",
    published_strain: str = "IMK1010",
    source_paper: str | None = None,
) -> pd.DataFrame:
    """Join predictions to published growth/no-growth and flag disagreements.

    Rows where the published phenotype is blank are reported as "unknown"
    rather than counted as agreement -- an unfilled cell is not evidence.

    Args:
        source_paper: restrict the comparison to one paper. This matters:
            de Valk et al. (2022) report growth of IMK1010 on 2% glycerol
            after 8 days, while Nummela et al. (2026) report no growth of the
            empty-vector control on 1% glycerol after 6 days. The two are not
            necessarily contradictory -- a slow passive route can be visible
            at one concentration and incubation time and not another -- but
            FBA cannot distinguish them, so we never average the two. We score
            against each source separately and report both.
    """
    published = phenotypes[phenotypes["strain"] == published_strain]
    if source_paper is not None:
        published = published[published["source_paper"] == source_paper]
    lookup = dict(zip(published["substrate"], published["published_grows"]))

    out = panel.copy()
    out["published"] = [str(lookup.get(s, "") or "").strip() for s in out["substrate"]]
    predicted = out[f"{strain} grows"]

    verdicts = []
    for pred, obs in zip(predicted, out["published"]):
        if pred == "not representable":
            # Yeast9 has no uptake route for this substrate under any strain,
            # so it can only ever predict "no growth". Scoring that as
            # agreement with an observed no-growth would be vacuous, and as
            # disagreement would be unfair. It is its own category.
            verdicts.append("not representable")
        elif obs not in {"yes", "no"}:
            verdicts.append("unknown")
        elif pred == obs:
            verdicts.append("agrees")
        else:
            verdicts.append("DISAGREES")
    out["agreement"] = verdicts
    return out


def uptake_survivor_table(
    wild_type: cobra.Model,
    deleted_orfs: list[str],
    substrates: dict[str, str] | None = None,
    orf_to_name: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Which genes keep each substrate's uptake open after the deletions.

    This is the evidence table for the project's main finding, so it is a
    result in its own right rather than a debugging aid.
    """
    substrates = substrates or DEVALK_SUBSTRATES
    orf_to_name = orf_to_name or {}

    records: list[dict[str, str]] = []
    for published_name, model_name in substrates.items():
        if not model_name:
            records.append({"substrate": published_name, "carbon_source": "",
                            "uptake_reaction": "", "surviving_genes": "",
                            "uptake_blocked_in_IMK1010": "n/a - absent from Yeast9"})
            continue

        # A condition may offer several carbon sources (ethanol + glycerol);
        # report the routes for each one separately.
        for single in [model_name] if isinstance(model_name, str) else model_name:
            survivors = strains.surviving_uptake_genes(wild_type, single, deleted_orfs)
            if not survivors:
                records.append({"substrate": published_name, "carbon_source": single,
                                "uptake_reaction": "", "surviving_genes": "",
                                "uptake_blocked_in_IMK1010": "yes"})
                continue
            for reaction_id, orfs in survivors.items():
                records.append({
                    "substrate": published_name,
                    "carbon_source": single,
                    "uptake_reaction": reaction_id,
                    "surviving_genes": ", ".join(orf_to_name.get(o, o) for o in orfs)
                    or "(no gene rule: passive diffusion)",
                    "uptake_blocked_in_IMK1010": "no",
                })

    return pd.DataFrame.from_records(records)


def write_csv(frame: pd.DataFrame, name: str, results_dir: pathlib.Path = RESULTS_DIR) -> pathlib.Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / name
    frame.to_csv(path, index=False)
    return path
