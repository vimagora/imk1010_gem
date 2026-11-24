"""Defining growth media in Yeast9.

The point of this module is to make one statement true and checkable:
*the model has exactly one carbon source open, and it is the one we asked for.*

Yeast9 ships with a sensible aerobic minimal medium already in place: ammonium
as nitrogen source, phosphate, sulphate, water, protons, oxygen and a set of
trace metals are open for uptake, and the only carbon source open is glucose.
So we do not build a medium from scratch; we take Yeast9's medium as given and
swap the carbon source. That keeps us honest (we are not inventing a medium
composition) and keeps the diff between conditions down to one reaction.

Conventions in Yeast9 that this module relies on, all verified against
v9.1.1 rather than assumed:

* An exchange reaction has a single metabolite with stoichiometry -1, e.g.
  ``r_1714: s_0565 -->`` is the D-glucose exchange.
* Flux is negative for uptake and positive for secretion. So opening an
  uptake means setting a *negative* lower bound, and closing one means
  setting the lower bound to 0 while leaving the upper bound alone, which
  still permits secretion.
"""

from __future__ import annotations

import cobra

# Yeast9's own default bound on glucose uptake (mmol gDW^-1 h^-1). We reuse it
# rather than pick a number of our own, so the default is traceable to the
# model and not to us. It is deliberately conservative: growth/no-growth
# conclusions do not depend on it, only predicted rates do.
DEFAULT_UPTAKE = 1.0


def exchange_reactions(model: cobra.Model) -> list[cobra.Reaction]:
    """Return every exchange reaction in the model.

    We identify them structurally (one metabolite, coefficient -1) rather than
    by an ID prefix or name, because reaction naming is not a contract.
    """
    found = []
    for reaction in model.reactions:
        if len(reaction.metabolites) == 1:
            (coefficient,) = reaction.metabolites.values()
            if coefficient == -1:
                found.append(reaction)
    return found


def exchanged_metabolite(reaction: cobra.Reaction) -> cobra.Metabolite:
    """The single metabolite an exchange reaction acts on."""
    (metabolite,) = reaction.metabolites
    return metabolite


def carbon_exchanges(model: cobra.Model) -> list[cobra.Reaction]:
    """Exchange reactions whose metabolite contains carbon.

    Determined from each metabolite's chemical formula, so the list cannot
    drift out of step with the model. Note this includes CO2 and fermentation
    products; that is intended, because we want to close every route by which
    carbon could *enter* the cell.
    """
    carbon = []
    for reaction in exchange_reactions(model):
        elements = exchanged_metabolite(reaction).elements or {}
        if elements.get("C", 0) > 0:
            carbon.append(reaction)
    return carbon


def find_exchange(model: cobra.Model, metabolite_name: str) -> cobra.Reaction:
    """Find the exchange reaction for a metabolite given its *name*.

    Yeast9 metabolite IDs are opaque (``s_0565``), so we look substrates up by
    name and resolve the ID here, once. Matching is case-insensitive and exact
    on the name, and we fail loudly on anything ambiguous rather than silently
    picking one -- a wrong exchange would quietly invalidate every result.
    """
    wanted = metabolite_name.strip().lower()
    matches = [
        reaction
        for reaction in exchange_reactions(model)
        if exchanged_metabolite(reaction).name.strip().lower() == wanted
    ]
    if not matches:
        raise KeyError(
            f"no exchange reaction for a metabolite named {metabolite_name!r}. "
            "Check the name against the model: Yeast9 uses some synonyms "
            "(D-sorbitol is stored as 'D-glucitol')."
        )
    if len(matches) > 1:
        raise KeyError(
            f"{metabolite_name!r} matches several exchanges: "
            f"{[r.id for r in matches]}. Disambiguate by reaction ID."
        )
    return matches[0]


def close_all_carbon_sources(model: cobra.Model) -> None:
    """Block uptake of every carbon-containing metabolite.

    Secretion is left untouched: the cell must still be able to excrete CO2,
    ethanol and so on, or nothing will grow.
    """
    for reaction in carbon_exchanges(model):
        reaction.lower_bound = 0.0


def set_carbon_source(
    model: cobra.Model,
    metabolite_name: str,
    uptake: float = DEFAULT_UPTAKE,
    per_carbon: bool = False,
) -> cobra.Reaction:
    """Make ``metabolite_name`` the model's only carbon source.

    Modifies ``model`` in place, so call it inside ``with model:`` to have
    cobrapy revert the change afterwards.

    Args:
        model: the model to modify.
        metabolite_name: name of the extracellular substrate, e.g. "glycerol".
        uptake: maximum uptake rate, mmol gDW^-1 h^-1 (a positive number).
        per_carbon: if True, treat ``uptake`` as a budget of carbon *atoms*
            rather than of molecules, so that substrates of different size are
            compared on equal carbon. Glycerol (C3) then gets twice the molar
            allowance of glucose (C6). Off by default because growth/no-growth
            results do not need it, but worth switching on before comparing
            predicted rates between substrates.

    Returns:
        The exchange reaction that was opened.
    """
    if uptake < 0:
        raise ValueError(f"uptake must be a positive magnitude, got {uptake}")

    exchange = find_exchange(model, metabolite_name)

    close_all_carbon_sources(model)

    if per_carbon:
        n_carbon = exchanged_metabolite(exchange).elements.get("C", 0)
        if not n_carbon:
            raise ValueError(f"{metabolite_name!r} contains no carbon")
        uptake = uptake / n_carbon

    # Negative lower bound = uptake allowed, up to this rate.
    exchange.lower_bound = -abs(uptake)
    return exchange


def open_carbon_sources(model: cobra.Model) -> list[str]:
    """Names of the carbon sources currently open for uptake.

    Used by the tests to assert the medium is what we think it is.
    """
    return sorted(
        exchanged_metabolite(reaction).name
        for reaction in carbon_exchanges(model)
        if reaction.lower_bound < 0
    )
