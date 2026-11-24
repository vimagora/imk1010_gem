"""Figures for Parts A, B and C.

Both figures are categorical state matrices, not magnitude heatmaps: each cell
holds one of four *states* (growth, no growth, not reported, not representable),
so each state gets one fill rather than a value ramp.

Three channels carry every state, so none of them is load-bearing alone:

* fill        -- dark blue for growth, near-surface tints for the rest
* glyph       -- printed in every cell (+ / - / ? / 0)
* texture     -- 45 degree hatch for "not representable", 135 degree for a
                 disagreement, following the texture rule of 45/135 only

That matters here because these figures are static images destined for a README
and possibly a printed PDF: there is no hover layer and no table view to fall
back on, so the state is printed in the cell. It also means the figures stay
readable in greyscale and under colour-vision deficiency, and nothing depends on
a red/green distinction.

Every number is read from results/*.csv. Run the analyses first:

    python scripts/run_part_a.py
    python scripts/run_part_b.py
    python scripts/run_part_c.py
    python scripts/make_figures.py
"""

from __future__ import annotations

import pathlib

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
FIGURES = RESULTS / "figures"

# --- palette (validated with the dataviz skill's validate_palette.js) --------
# "#2a78d6,#d03b3b" light mode: all checks pass, worst adjacent CVD dE 23.8
# (protan), normal-vision dE 31.6, both >= 3:1 against the light surface.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
HAIRLINE = "#c3c2b7"
GROWTH = "#2a78d6"
NO_GROWTH = "#eeede8"
ABSENT = "#dedcd3"
CRITICAL = "#d03b3b"

# state -> (fill, glyph, glyph colour, hatch)
STATES = {
    "yes": (GROWTH, "+", "#ffffff", None),
    "no": (NO_GROWTH, "−", INK_SECONDARY, None),
    "not reported": (SURFACE, "?", INK_MUTED, None),
    "not representable": (ABSENT, "⌀", INK_SECONDARY, "////"),
}

# Figure 3 also has to show which *components* a construct carries, which is a
# different question from whether it grew. Those cells therefore use a separate
# neutral vocabulary: dark grey for a component that is present, and glyphs
# that do not collide with the growth glyphs (+ / −). Growth blue is reserved
# for growth alone, so the two column groups can never be misread as each other.
COMPONENT = "#55544f"
COMPONENT_STATES = {
    "intracellular": (COMPONENT, "IN", "#ffffff", None),
    "secreted": (COMPONENT, "OUT", "#ffffff", "////"),
    "present": (COMPONENT, "●", "#ffffff", None),
    "absent": (NO_GROWTH, "○", INK_MUTED, None),
}

# A fixed, meaningful row order: hexoses, disaccharides, polyols, then the
# combined plate. Rows absent from a panel are dropped, never reordered.
ROW_ORDER = [
    "D-glucose", "D-fructose", "D-mannose", "D-galactose",
    "sucrose", "maltose",
    "D-sorbitol", "D-xylitol", "D-mannitol",
    "glycerol", "ethanol + glycerol",
]

CELL = 0.86  # cell size, leaving a surface gap between cells


def _order(substrates: list[str]) -> list[str]:
    known = [s for s in ROW_ORDER if s in substrates]
    return known + sorted(set(substrates) - set(known))


def _state_of(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"yes", "no", "not reported", "not representable"}:
        return text
    return "not reported"


def draw_cell(ax, x: int, y: int, state: str, disagrees: bool = False,
              states: dict | None = None) -> None:
    """Draw one state cell, with a disagreement overlay if needed.

    ``states`` selects the vocabulary: the growth states by default, or
    COMPONENT_STATES for Figure 3's construct columns.
    """
    fill, glyph, glyph_colour, hatch = (states or STATES)[state]
    glyph_size = 11 if len(glyph) == 1 else 7.5

    ax.add_patch(
        FancyBboxPatch(
            (x - CELL / 2, y - CELL / 2), CELL, CELL,
            boxstyle="round,pad=0,rounding_size=0.1",
            facecolor=fill,
            edgecolor=HAIRLINE if state != "not reported" else INK_MUTED,
            linewidth=0.6,
            linestyle="--" if state == "not reported" else "-",
            hatch=hatch,
            zorder=2,
        )
    )
    if hatch:
        # Hatch colour follows the edge colour in matplotlib, so redraw the
        # outline on top at the intended hairline weight.
        ax.add_patch(
            FancyBboxPatch(
                (x - CELL / 2, y - CELL / 2), CELL, CELL,
                boxstyle="round,pad=0,rounding_size=0.1",
                facecolor="none", edgecolor=HAIRLINE, linewidth=0.6, zorder=3,
            )
        )

    glyph_text = ax.text(x, y, glyph, ha="center", va="center", color=glyph_colour,
                         fontsize=glyph_size, fontweight="bold", zorder=5)
    if hatch:
        # Hatch lines cross the glyph and can turn "OUT" into "ØUT". A stroke
        # in the cell's own fill colour keeps the glyph legible.
        glyph_text.set_path_effects(
            [patheffects.withStroke(linewidth=2.6, foreground=fill)]
        )

    if disagrees:
        ax.add_patch(
            FancyBboxPatch(
                (x - CELL / 2, y - CELL / 2), CELL, CELL,
                boxstyle="round,pad=0,rounding_size=0.1",
                facecolor="none", edgecolor=CRITICAL, linewidth=2.2,
                hatch="\\\\\\\\", zorder=4,
            )
        )


def _frame(ax, columns: list[str], rows: list[str], title: str, subtitle: str,
           header_offset: int = 34,
           group_headers: list[tuple[str, float, float]] | None = None,
           group_y: float = -1.12) -> None:
    """Set up one matrix panel.

    Titles are positioned in *offset points* from the top of the axes rather
    than in axes fractions, so panels with different row counts keep identical
    spacing instead of drifting into their column headers.
    """
    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)  # first row at the top
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(columns, fontsize=8, color=INK_SECONDARY, linespacing=1.5)
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=9, color=INK)
    ax.tick_params(length=0, pad=5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor(SURFACE)

    page_left = ("figure fraction", "axes fraction")
    ax.annotate(title, xy=(0.03, 1), xycoords=page_left,
                xytext=(0, header_offset + 16), textcoords="offset points",
                fontsize=10.5, fontweight="bold", color=INK,
                ha="left", va="bottom", annotation_clip=False)
    ax.annotate(subtitle, xy=(0.03, 1), xycoords=page_left,
                xytext=(0, header_offset + 2), textcoords="offset points",
                fontsize=8, color=INK_SECONDARY, linespacing=1.5,
                ha="left", va="bottom", annotation_clip=False)

    if group_headers:
        for label, start, end in group_headers:
            ax.annotate(
                label, xy=((start + end) / 2, group_y), xycoords="data",
                xytext=(0, 6), textcoords="offset points",
                ha="center", va="bottom", fontsize=9.5, fontweight="bold",
                color=INK, annotation_clip=False,
            )
            ax.plot([start - 0.38, end + 0.38], [group_y, group_y],
                    color=HAIRLINE, linewidth=1.0, clip_on=False,
                    solid_capstyle="butt")


def legend_handles(include_disagreement: bool = True) -> list:
    labels = {
        "yes": "growth",
        "no": "no growth",
        "not reported": "not reported in the paper",
        "not representable": "not representable in Yeast9",
    }
    handles = []
    for state, label in labels.items():
        fill, glyph, _, hatch = STATES[state]
        handles.append(
            Line2D([], [], marker="s", linestyle="none", markersize=11,
                   markerfacecolor=fill, markeredgecolor=HAIRLINE,
                   label=f"{glyph}   {label}")
        )
    if include_disagreement:
        handles.append(
            Line2D([], [], marker="s", linestyle="none", markersize=11,
                   markerfacecolor="none", markeredgecolor=CRITICAL,
                   markeredgewidth=2.2,
                   label="prediction disagrees with the published phenotype")
        )
    return handles


# --------------------------------------------------------------------------
# Figure 1: Part A
# --------------------------------------------------------------------------

def figure_part_a(panel: pd.DataFrame) -> pathlib.Path:
    strains = ["wild type", "IMK1010", "IMK1010 repaired"]
    sources = [
        ("de Valk et al. 2022", "deValk2022", "IMK1010 (Sugar0) on SM plates"),
        ("Nummela et al. 2026", "Nummela2026", "IMK1010 + empty vector, solid medium"),
    ]

    # Row counts differ between panels, so size each panel by its own rows and
    # keep the cells square.
    subsets = {}
    for _, paper, _ in sources:
        rows = panel[panel[f"[{paper}] published"].isin(["yes", "no"])]
        subsets[paper] = rows.set_index("substrate").loc[
            _order(list(rows["substrate"]))
        ]

    heights = [len(subsets[p]) for _, p, _ in sources]
    fig = plt.figure(figsize=(7.6, 8.8), facecolor=SURFACE)
    grid = fig.add_gridspec(
        2, 1, height_ratios=heights, hspace=0.62,
        left=0.24, right=0.94, top=0.80, bottom=0.13,
    )

    columns = ["published"] + strains

    for index, (label, paper, detail) in enumerate(sources):
        ax = fig.add_subplot(grid[index])
        table = subsets[paper]
        rows = list(table.index)

        # The agreement score belongs under the column it describes, not in a
        # subtitle that would run off the page.
        scores = {}
        for strain in ["IMK1010", "IMK1010 repaired"]:
            verdicts = table[f"[{paper}] {strain}"]
            scoreable = verdicts[verdicts != "not representable"]
            scores[strain] = f"{(scoreable == 'agrees').sum()}/{len(scoreable)}"

        column_labels = [
            "published",
            "wild type\n(context)",
            f"IMK1010\n(faithful)\n{scores['IMK1010']}",
            f"IMK1010\n(repaired)\n{scores['IMK1010 repaired']}",
        ]

        _frame(
            ax, column_labels, rows,
            f"{'AB'[index]}.  Measured against {label}", detail,
            header_offset=40,
        )

        for row_index, substrate in enumerate(rows):
            for col_index, column in enumerate(columns):
                if column == "published":
                    state = _state_of(table.loc[substrate, f"[{paper}] published"])
                    disagrees = False
                else:
                    state = _state_of(table.loc[substrate, f"{column} grows"])
                    verdict = str(table.loc[substrate, f"[{paper}] {column}"]) \
                        if column != "wild type" else ""
                    disagrees = verdict == "DISAGREES"
                draw_cell(ax, col_index, row_index, state, disagrees)

    fig.text(0.03, 0.975, "Figure 1  Predicted vs published growth of IMK1010 in Yeast9",
             fontsize=12, fontweight="bold", color=INK, va="top", ha="left")
    fig.text(
        0.03, 0.945,
        "The repaired strain additionally deletes the genes Yeast9 retains in sugar uptake gene rules\n"
        "(HXT17, CIN10, GUP1/GUP2, FPS1), derived from the model's own rules. It exists in no laboratory.",
        fontsize=8.5, color=INK_SECONDARY, va="top", ha="left", linespacing=1.5,
    )
    fig.legend(
        handles=legend_handles(), loc="lower left", bbox_to_anchor=(0.03, 0.012),
        frameon=False, fontsize=8.5, ncol=2, handletextpad=0.6,
        columnspacing=2.0, labelcolor=INK_SECONDARY,
    )

    return _save(fig, "figure1_part_a_carbon_panel")


# --------------------------------------------------------------------------
# Figure 2: Part B
# --------------------------------------------------------------------------

def figure_part_b(table: pd.DataFrame) -> pathlib.Path:
    porters = sorted(table["porter"].unique())
    rows = _order(list(table["substrate"].unique()))

    fig = plt.figure(figsize=(7.6, 6.4), facecolor=SURFACE)
    ax = fig.add_subplot(
        fig.add_gridspec(1, 1, left=0.34, right=0.95, top=0.66, bottom=0.18)[0]
    )

    columns, labels, groups = [], [], []
    for index, porter in enumerate(porters):
        start = index * 2
        columns += [(porter, "published"), (porter, "predicted")]
        labels += ["published", "predicted"]
        groups.append((porter, float(start), float(start + 1)))

    _frame(ax, labels, rows,
           "Published substrate specificity vs model prediction",
           "Background: IMK1010 repaired, which grows on none of these substrates on its own",
           header_offset=54, group_headers=groups)

    lookup = table.set_index(["porter", "substrate"])
    for row_index, substrate in enumerate(rows):
        for col_index, (porter, kind) in enumerate(columns):
            record = lookup.loc[(porter, substrate)]
            state = _state_of(record["published" if kind == "published" else "predicted"])
            disagrees = str(record["verdict"]) == "DISAGREES"
            draw_cell(ax, col_index, row_index, state,
                      disagrees and kind == "predicted")

    # Make the discriminating result explicit: only GlpC carries the hexoses.
    hexoses = [r for r in ["D-glucose", "D-fructose", "D-mannose"] if r in rows]
    if hexoses:
        top = rows.index(hexoses[0]) - 0.5
        bottom = rows.index(hexoses[-1]) + 0.5
        ax.plot([-1.78, -1.78], [top, bottom], color=CRITICAL, linewidth=2.2,
                clip_on=False, solid_capstyle="butt")
        ax.annotate(
            "hexoses:\nGlpC only", xy=(-1.88, (top + bottom) / 2),
            ha="right", va="center", fontsize=8.5, color=CRITICAL,
            fontweight="bold", annotation_clip=False, linespacing=1.5,
        )

    claimed = table[table["published"] == "yes"]
    scoreable = claimed[claimed["predicted"] != "not representable"]
    n_agree = (scoreable["verdict"] == "agrees").sum()

    fig.text(0.03, 0.975,
             "Figure 2  Aspergillus niger Glp porters added to IMK1010 in silico",
             fontsize=12, fontweight="bold", color=INK, va="top", ha="left")
    fig.text(
        0.03, 0.945,
        f"Of {len(claimed)} published porter-substrate pairs, {len(scoreable)} are representable in Yeast9 "
        f"and {n_agree} are reproduced.\n"
        "The informative half is negative: GlpA and GlpB do not restore hexose growth, GlpC does.\n"
        "\"Not reported\" is not evidence of absence — the paper reports no result for those pairs.",
        fontsize=8.5, color=INK_SECONDARY, va="top", ha="left",
    )
    fig.legend(
        handles=legend_handles(include_disagreement=False), loc="lower left",
        bbox_to_anchor=(0.03, 0.012), frameon=False, fontsize=8.5, ncol=2,
        handletextpad=0.6, columnspacing=2.0, labelcolor=INK_SECONDARY,
    )

    return _save(fig, "figure2_part_b_glp_porters")


# --------------------------------------------------------------------------
# Figure 3: Part C
# --------------------------------------------------------------------------

def _components(construct: str) -> tuple[str, str, str]:
    """Read a construct description into (beta-glucosidase, transporter, detail).

    Parsed from the construct strings written by run_part_c.py rather than
    hardcoded, so the figure cannot drift out of step with the analysis.
    """
    text = construct.lower()
    if "secreted bgl" in text:
        bgl = "secreted"
    elif "intracellular bgl" in text:
        bgl = "intracellular"
    else:
        bgl = "absent"

    if "cta" in text:
        transporter = "present"
        detail = "CtA, H⁺ symport" if "symport" in text else "CtA, facilitated"
    else:
        transporter = "absent"
        detail = "—"
    return bgl, transporter, detail


def figure_part_c(table: pd.DataFrame) -> pathlib.Path:
    # Keep the order run_part_c.py produced: it is the argument of the
    # experiment, from "nothing added" through to the two controls.
    rows = list(table.index)
    # Short column names: "enzyme" rather than the full name, which would not
    # fit a cell width. The legend spells out what the enzyme is.
    labels = ["enzyme", "transporter", "expected", "observed"]

    # Row labels carry the transporter variant, so the matrix needs no extra
    # text column and the right margin is free for the rate and the findings.
    row_labels = []
    for key in rows:
        record = table.loc[key]
        _, _, detail = _components(str(record["construct"]))
        suffix = f"\n({detail.split(', ')[-1]})" if detail != "—" else ""
        row_labels.append(f"{record['strain']}{suffix}")

    fig = plt.figure(figsize=(8.0, 6.6), facecolor=SURFACE)
    ax = fig.add_subplot(
        fig.add_gridspec(1, 1, left=0.31, right=0.68, top=0.62, bottom=0.27)[0]
    )

    _frame(ax, labels, row_labels,
           "Cellobiose needs both a transporter and a hydrolase",
           "Which one is missing is the experiment; the strain supplies only the enzyme",
           header_offset=52,
           group_headers=[("construct", 0.0, 1.0),
                          ("growth on cellobiose", 2.0, 3.0)],
           group_y=-1.62)
    # Two-line row labels need a smaller size than the single-line labels the
    # other figures use, or consecutive rows touch.
    ax.set_yticklabels(row_labels, fontsize=8.3, color=INK, linespacing=1.4)

    markers = {2: "1", 3: "1", 4: "1", 5: "2", 6: "2"}

    for row_index, key in enumerate(rows):
        record = table.loc[key]
        bgl, transporter, _ = _components(str(record["construct"]))

        draw_cell(ax, 0, row_index, bgl, states=COMPONENT_STATES)
        draw_cell(ax, 1, row_index, transporter, states=COMPONENT_STATES)
        draw_cell(ax, 2, row_index, _state_of(record["expected"]))
        draw_cell(ax, 3, row_index, _state_of(record["grows"]),
                  disagrees=str(record["matches_expectation"]).upper() == "NO")

        # The symport/facilitated difference is invisible as a state and only
        # shows up as a rate, so the rate is printed alongside.
        rate = float(record["growth_rate"])
        ax.annotate(f"{rate:.4f}" if rate else "—",
                    xy=(4.45, row_index), xycoords="data",
                    ha="right", va="center", fontsize=8.5,
                    color=INK if rate else INK_MUTED,
                    fontweight="bold" if rate else "normal",
                    annotation_clip=False)
        if row_index in markers:
            ax.annotate(markers[row_index], xy=(4.95, row_index),
                        xycoords="data", ha="center", va="center",
                        fontsize=8, color=CRITICAL, fontweight="bold",
                        annotation_clip=False)

    ax.set_xlim(-0.5, 3.5)
    ax.annotate("growth (h⁻¹)", xy=(4.45, -1.02), xycoords="data", ha="right",
                va="center", fontsize=7.5, color=INK_SECONDARY,
                annotation_clip=False)

    # Brackets tying each numbered note to the rows that carry it.
    for top, bottom, label in [(2, 4, "1"), (5, 6, "2")]:
        ax.plot([4.72, 4.72], [top - 0.44, bottom + 0.44], color=CRITICAL,
                linewidth=2.0, clip_on=False, solid_capstyle="butt")

    matches = (table["matches_expectation"].str.lower() == "yes").sum()

    fig.text(0.03, 0.975,
             "Figure 3  A cellobiose platform strain: IMK1010 + NcBGL",
             fontsize=12, fontweight="bold", color=INK, va="top", ha="left")
    fig.text(
        0.03, 0.945,
        "S. cerevisiae lacks both a cellobiose transporter and a beta-glucosidase, so two things are missing.\n"
        f"Supplying only the enzyme isolates transport as the variable. {matches}/{len(table)} results match expectation.\n"
        "Cellobiose is absent from Yeast9 entirely: this part extends the model, and is reported as an extension.",
        fontsize=8.5, color=INK_SECONDARY, va="top", ha="left", linespacing=1.5,
    )

    # The two findings, full width so they are not cramped into a margin.
    for offset, (number, note) in enumerate([
        ("1", "The enzyme alone does not grow; adding CtA does. "
              "Transport is the limiting step, not the hydrolase."),
        ("2", "The same enzyme secreted instead of intracellular rescues wild type "
              "but not IMK1010:\n     the glucose it releases stays outside, and "
              "IMK1010 has no hexose transporter left."),
    ]):
        fig.text(0.03, 0.225 - offset * 0.055, number, fontsize=8.5,
                 color=CRITICAL, fontweight="bold", va="top", ha="left")
        fig.text(0.055, 0.225 - offset * 0.055, note, fontsize=8.5,
                 color=INK_SECONDARY, va="top", ha="left", linespacing=1.5)

    handles = [
        Line2D([], [], marker="s", linestyle="none", markersize=11,
               markerfacecolor=COMPONENT_STATES[state][0],
               markeredgecolor=HAIRLINE,
               label=f"{COMPONENT_STATES[state][1]}   {label}")
        for state, label in [("intracellular", "enzyme, intracellular"),
                             ("secreted", "enzyme, secreted"),
                             ("present", "component present"),
                             ("absent", "component absent")]
    ] + [
        Line2D([], [], marker="s", linestyle="none", markersize=11,
               markerfacecolor=GROWTH, markeredgecolor=HAIRLINE, label="+   growth"),
        Line2D([], [], marker="s", linestyle="none", markersize=11,
               markerfacecolor=NO_GROWTH, markeredgecolor=HAIRLINE,
               label="−   no growth"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.03, 0.005),
               frameon=False, fontsize=8.5, ncol=3, handletextpad=0.6,
               columnspacing=1.8, labelcolor=INK_SECONDARY)

    return _save(fig, "figure3_part_c_cellobiose")


def _save(fig, stem: str) -> pathlib.Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    png = FIGURES / f"{stem}.png"
    fig.savefig(png, dpi=200, facecolor=SURFACE)
    # Suppress the PDF creation timestamp: without this the file differs on
    # every run and shows up as a spurious change in version control.
    fig.savefig(FIGURES / f"{stem}.pdf", facecolor=SURFACE,
                metadata={"CreationDate": None})  # vector, for print
    plt.close(fig)
    return png


def main() -> None:
    panel = pd.read_csv(RESULTS / "carbon_panel.csv")
    glp = pd.read_csv(RESULTS / "glp_predictions.csv")
    cellobiose = pd.read_csv(RESULTS / "cellobiose_panel.csv")

    for path in [figure_part_a(panel), figure_part_b(glp),
                 figure_part_c(cellobiose)]:
        print(f"wrote {path}")
        print(f"wrote {path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
