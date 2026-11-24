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
# Where the two differ, the published name is kept for reporting.
DEVALK_SUBSTRATES = {
    "D-glucose": "D-glucose",
    "D-fructose": "D-fructose",
    "D-galactose": "D-galactose",
    "sucrose": "sucrose",
    "maltose": "maltose",
    "ethanol": "ethanol",
    "glycerol": "glycerol",
}


def growth_on(
    model: cobra.Model,
    substrate: str,
    uptake: float = media.DEFAULT_UPTAKE,
    per_carbon: bool = False,
) -> float:
    """Predicted growth rate (h^-1) with ``substrate`` as the only carbon source.

    Returns 0.0 for an infeasible problem, which is what "cannot grow" looks
    like to the solver.
    """
    with model:
        media.set_carbon_source(model, substrate, uptake=uptake, per_carbon=per_carbon)
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
        row: dict[str, object] = {
            "substrate": published_name,
            "model_metabolite": model_name,
            "exchange": media.find_exchange(wild_type, model_name).id,
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
) -> pd.DataFrame:
    """Join predictions to published growth/no-growth and flag disagreements.

    Rows where the published phenotype is blank are reported as "unknown"
    rather than counted as agreement -- an unfilled cell is not evidence.
    """
    published = phenotypes[phenotypes["strain"] == published_strain]
    lookup = dict(zip(published["substrate"], published["published_grows"]))

    out = panel.copy()
    out["published"] = [str(lookup.get(s, "") or "").strip() for s in out["substrate"]]
    predicted = out[f"{strain} grows"]

    verdicts = []
    for pred, obs in zip(predicted, out["published"]):
        if obs not in {"yes", "no"}:
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

    records = []
    for published_name, model_name in substrates.items():
        survivors = strains.surviving_uptake_genes(wild_type, model_name, deleted_orfs)
        if not survivors:
            records.append(
                {
                    "substrate": published_name,
                    "uptake_reaction": "",
                    "surviving_genes": "",
                    "uptake_blocked_in_IMK1010": "yes",
                }
            )
            continue
        for reaction_id, orfs in survivors.items():
            records.append(
                {
                    "substrate": published_name,
                    "uptake_reaction": reaction_id,
                    "surviving_genes": ", ".join(orf_to_name.get(o, o) for o in orfs)
                    or "(no gene rule: passive diffusion)",
                    "uptake_blocked_in_IMK1010": "no",
                }
            )
    return pd.DataFrame.from_records(records)


def write_csv(frame: pd.DataFrame, name: str, results_dir: pathlib.Path = RESULTS_DIR) -> pathlib.Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / name
    frame.to_csv(path, index=False)
    return path
