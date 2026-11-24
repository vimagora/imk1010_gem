"""Loading Yeast9 and recovering the gene names it hides.

Yeast9 stores each gene's systematic ORF name as the SBML gene product ID
(``YHR094C``) and its common name in the ``fbc:label`` attribute (``HXT1``).
COBRApy does not expose that label, so ``gene.name`` comes back empty and
there is no way to look a gene up by its common name using cobrapy alone.

Since our gene list from the papers is written in common names, we parse the
labels straight out of the SBML. The mapping therefore comes from the model
file itself rather than from anybody's memory -- which is the point.
"""

from __future__ import annotations

import pathlib
import re
import warnings

import cobra

MODEL_PATH = pathlib.Path(__file__).resolve().parents[2] / "data" / "model" / "yeast-GEM.xml"

# <fbc:geneProduct ... fbc:id="YHR094C" fbc:label="HXT1">
_GENE_PRODUCT = re.compile(r'<fbc:geneProduct[^>]*fbc:id="([^"]+)"[^>]*fbc:label="([^"]+)"')


def load_model(path: pathlib.Path | str = MODEL_PATH) -> cobra.Model:
    """Read the Yeast9 SBML file.

    SBML parsing emits warnings about unsupported annotation packages that are
    irrelevant here, so they are suppressed to keep output readable.
    """
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run: python scripts/download_model.py"
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cobra.io.read_sbml_model(str(path))


def gene_name_map(path: pathlib.Path | str = MODEL_PATH) -> dict[str, str]:
    """Map systematic ORF ID -> common gene name, read from the SBML labels."""
    text = pathlib.Path(path).read_text(encoding="utf-8")
    return dict(_GENE_PRODUCT.findall(text))


def common_name_to_orf(path: pathlib.Path | str = MODEL_PATH) -> dict[str, str]:
    """Map upper-cased common gene name -> systematic ORF ID."""
    mapping: dict[str, str] = {}
    for orf, label in gene_name_map(path).items():
        # First label wins; duplicates would be a model-level oddity worth seeing.
        mapping.setdefault(label.upper(), orf)
    return mapping


def describe_gpr(rule: str, orf_to_name: dict[str, str]) -> str:
    """Rewrite a gene-protein-reaction rule into common gene names.

    Gene rules in Yeast9 read ``YHR094C or YMR011W or ...``, which is unreadable
    when the question is "is HXT17 in this rule?". This makes them legible for
    printing and for the README.
    """
    if not rule:
        return "(no gene rule)"
    return re.sub(
        r"\b(Y[A-Z]{2}\d{3}[WC](?:-[A-Z])?|Q\d{4})\b",
        lambda match: orf_to_name.get(match.group(0), match.group(0)),
        rule,
    )
