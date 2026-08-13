"""Mapping our gene ids onto the Wang et al. DE gene lists, for external validation.

Three naming systems are involved:

* ours          ``Glyma.01G161700.Wm82.a4.v1``   (Wm82.a4.v1, with an assembly suffix)
* correspondence``Glyma.01g000100`` -> ``Glyma01g00210``  (Wm82.a4.v1 -> Glyma 1.1)
* Wang          ``GLYMA11G21190``, sometimes ``GLYMA19G01050.13`` (Glyma 1.1, uppercased,
                occasionally carrying an isoform suffix)

Everything is compared uppercased and suffix-stripped.

Why this matters for the question at hand: the Wang lists are an *independent* bulk
RNA-seq measurement of the 24h heat, drought and combined responses. If a processing
choice were genuinely destroying real Heat/Drought signal, the strategies that recover
more Wang genes are the ones preserving it. That turns "which pipeline is right" from
an argument about defaults into something with an external criterion.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import pandas as pd

from . import config as cfg

REFERENCE_DIR = cfg.REPO_ROOT / "Referenced_data"
CORRESPONDENCE = REFERENCE_DIR / "Wm82.a4.v1_to_Correspondence_Full.csv"
WANG_DIR = REFERENCE_DIR / "Wang_DE_gene_lists"

# Wang list -> the contrast in our analysis it corresponds to
WANG_SETS = {
    "Heat vs Control": "heat_24h_DE.csv",
    "Drought vs Control": "only_drought_24h_DE.csv",
    "HD vs Control": "drough_heat_24h_DE.csv",
}


def _norm(gene: str) -> str:
    """Uppercase, drop our assembly suffix and any trailing isoform number."""
    g = str(gene).strip().upper()
    g = g.replace(".WM82.A4.V1", "")
    g = re.sub(r"\.\d+$", "", g)          # GLYMA19G01050.13 -> GLYMA19G01050
    return g


@lru_cache(maxsize=1)
def a4_to_v1() -> dict[str, str]:
    """Wm82.a4.v1 gene -> Glyma 1.1 gene, normalised on both sides."""
    df = pd.read_csv(CORRESPONDENCE)
    df = df[df["Assembly"].astype(str).str.strip().str.lower() == "glyma 1.1"]
    out: dict[str, str] = {}
    for a4, v1 in zip(df["Wm82.a4.v1 Gene"], df["Corresponding Gene"]):
        out.setdefault(_norm(a4), _norm(v1))
    return out


def map_to_wang_space(genes) -> pd.Series:
    """Our var_names -> Glyma 1.1 ids (NaN where no correspondence exists)."""
    mapping = a4_to_v1()
    normed = [_norm(g) for g in genes]
    return pd.Series([mapping.get(g) for g in normed], index=list(genes), name="glyma_v1")


@lru_cache(maxsize=8)
def wang_set(contrast: str, min_probability: float = 0.0) -> frozenset[str]:
    path = WANG_DIR / WANG_SETS[contrast]
    df = pd.read_csv(path)
    if min_probability and "Probability" in df:
        df = df[df["Probability"].astype(float) >= min_probability]
    return frozenset(_norm(g) for g in df["GeneID"])


def wang_direction(contrast: str) -> dict[str, str]:
    df = pd.read_csv(WANG_DIR / WANG_SETS[contrast])
    col = [c for c in df.columns if c.startswith("Up-Down")][0]
    return {_norm(g): str(d) for g, d in zip(df["GeneID"], df[col])}


def enrichment(
    called: pd.Index,
    tested: pd.Index,
    contrast: str,
    min_probability: float = 0.0,
) -> dict:
    """Hypergeometric enrichment of a Wang set among our called genes.

    The background is restricted to genes that were actually *tested* and that have a
    correspondence entry -- using all genes would inflate the enrichment, since
    testability itself is not random.
    """
    from scipy import stats as sps

    ref = wang_set(contrast, min_probability)
    tested_map = map_to_wang_space(tested).dropna()
    called_map = map_to_wang_space(called).dropna()

    background = set(tested_map.values)
    ref_in_bg = ref & background
    hits = set(called_map.values) & ref_in_bg

    n_bg, n_ref, n_called, n_hit = len(background), len(ref_in_bg), len(set(called_map.values)), len(hits)
    if not n_bg or not n_ref or not n_called:
        return {"contrast": contrast, "n_called": n_called, "n_wang_in_background": n_ref,
                "n_recovered": n_hit, "recall": float("nan"), "fold_enrichment": float("nan"),
                "pvalue": float("nan")}

    expected = n_called * n_ref / n_bg
    pval = sps.hypergeom.sf(n_hit - 1, n_bg, n_ref, n_called)
    return {
        "contrast": contrast,
        "n_called": n_called,
        "n_wang_in_background": n_ref,
        "n_recovered": n_hit,
        "recall": n_hit / n_ref,
        "fold_enrichment": n_hit / expected if expected else float("nan"),
        "pvalue": float(pval),
    }
