"""Loading CellRanger / CellBender matrices and comparing two versions of one library."""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from . import config as cfg


# --------------------------------------------------------------------------- loading


def read_matrix(path: Path, var_names: str = "gene_ids") -> ad.AnnData:
    """Read a 10x mtx directory.

    ``gene_ids`` (not symbols) is the default because it is the only identifier
    guaranteed unique and stable across CellRanger versions -- symbol handling is
    exactly the kind of thing that changes between 9.0.1 and 10.1.0, and using it
    as the join key would manufacture differences that are not in the counts.
    """
    adata = sc.read_10x_mtx(path, var_names=var_names, make_unique=False)
    adata.var_names_make_unique()
    return adata


def read_metrics(outs: Path) -> pd.Series:
    df = pd.read_csv(outs / "metrics_summary.csv")
    row = df.iloc[0]
    return pd.Series({k: _to_number(v) for k, v in row.items()})


def _to_number(value):
    if isinstance(value, str):
        stripped = value.replace(",", "").strip()
        if stripped.endswith("%"):
            return float(stripped[:-1]) / 100.0
        try:
            return float(stripped)
        except ValueError:
            return value
    return value


def read_cellbender(h5: Path) -> ad.AnnData:
    """CellBender's filtered .h5. Falls back to scanpy's 10x reader."""
    adata = sc.read_10x_h5(h5)
    adata.var_names_make_unique()
    return adata


def annotate_organelles(adata: ad.AnnData) -> ad.AnnData:
    """Flag chloroplast / mitochondrial genes and add per-cell fractions.

    Uses ``var['gene_symbols']`` when the object was loaded by gene id, since the
    organelle contigs are identified by their symbol (``GlmaCp*`` / ``GlmaCt*``).
    """
    names = adata.var.get("gene_symbols", pd.Series(adata.var_names, index=adata.var_names))
    names = pd.Index(names.astype(str))
    adata.var["chloroplast"] = names.str.startswith(cfg.CHLOROPLAST_PREFIX)
    adata.var["mito"] = np.logical_or.reduce(
        [names.str.startswith(p) for p in cfg.MITO_PREFIXES]
    )
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["chloroplast", "mito"], inplace=True, log1p=False, percent_top=None
    )
    return adata


# --------------------------------------------------------------------------- comparing


def align(a: ad.AnnData, b: ad.AnnData) -> tuple[sparse.csr_matrix, sparse.csr_matrix, pd.Index, pd.Index]:
    """Restrict both matrices to their shared barcodes and genes, in the same order."""
    obs = a.obs_names.intersection(b.obs_names)
    var = a.var_names.intersection(b.var_names)
    A = a[obs, var].X
    B = b[obs, var].X
    return sparse.csr_matrix(A), sparse.csr_matrix(B), obs, var


def compare_matrices(a: ad.AnnData, b: ad.AnnData, label_a: str, label_b: str) -> dict:
    """Everything needed to say 'these are/aren't the same count matrix'."""
    shared_obs = a.obs_names.intersection(b.obs_names)
    shared_var = a.var_names.intersection(b.var_names)

    out = {
        f"n_cells_{label_a}": a.n_obs,
        f"n_cells_{label_b}": b.n_obs,
        f"n_genes_{label_a}": a.n_vars,
        f"n_genes_{label_b}": b.n_vars,
        "n_cells_shared": len(shared_obs),
        "n_genes_shared": len(shared_var),
        "cell_jaccard": len(shared_obs) / max(len(a.obs_names.union(b.obs_names)), 1),
        "gene_jaccard": len(shared_var) / max(len(a.var_names.union(b.var_names)), 1),
        f"total_umis_{label_a}": float(a.X.sum()),
        f"total_umis_{label_b}": float(b.X.sum()),
    }
    if not len(shared_obs) or not len(shared_var):
        return out

    A, B, _, _ = align(a, b)
    diff = (A - B).tocoo()
    nnz_union = int(((A != 0) + (B != 0)).nnz)

    out.update(
        {
            "identical": bool(diff.nnz == 0),
            "n_differing_entries": int(diff.nnz),
            "frac_entries_differing": diff.nnz / max(nnz_union, 1),
            "max_abs_diff": float(np.abs(diff.data).max()) if diff.nnz else 0.0,
            "sum_abs_diff": float(np.abs(diff.data).sum()),
            "shared_umis_" + label_a: float(A.sum()),
            "shared_umis_" + label_b: float(B.sum()),
            # On the shared cells only -- the whole-object totals above also fold in
            # cell-calling differences, which makes them useless for judging how much
            # a denoising step actually changed the counts.
            "frac_umis_removed_shared_cells": float(1 - B.sum() / A.sum()) if A.sum() else np.nan,
        }
    )

    # Direction of the disagreement. A pure read-depth difference (one of two
    # sequencing runs missing) makes ``b`` uniformly *lower*: essentially every
    # differing entry is down and no cell gains UMIs. A genuine alignment or
    # reference difference scatters in both directions. This is the column that
    # separates those two explanations.
    if diff.nnz:
        # ``diff`` is a - b, so a positive entry means b is the lower of the two.
        out["frac_diff_entries_lower_in_" + label_b] = float((diff.data > 0).mean())
    out["n_cells_higher_in_" + label_b] = int(
        (np.asarray(B.sum(axis=1)).ravel() > np.asarray(A.sum(axis=1)).ravel()).sum()
    )

    per_cell_a = np.asarray(A.sum(axis=1)).ravel()
    per_cell_b = np.asarray(B.sum(axis=1)).ravel()
    per_gene_a = np.asarray(A.sum(axis=0)).ravel()
    per_gene_b = np.asarray(B.sum(axis=0)).ravel()
    out["per_cell_umi_pearson"] = _corr(per_cell_a, per_cell_b)
    out["per_gene_umi_pearson"] = _corr(per_gene_a, per_gene_b)
    out["per_gene_umi_spearman"] = _corr(
        pd.Series(per_gene_a).rank().to_numpy(), pd.Series(per_gene_b).rank().to_numpy()
    )
    out["median_per_cell_ratio"] = float(
        np.median(per_cell_b[per_cell_a > 0] / per_cell_a[per_cell_a > 0])
    ) if (per_cell_a > 0).any() else np.nan
    return out


def _corr(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson r between two 1-D vectors.

    Written out element-wise rather than via ``np.corrcoef``, which routes
    through ``np.cov`` -> gemm. Same answer, but it keeps a whole-script
    dependency on a working BLAS out of what is only ever a correlation of two
    vectors -- on Windows an unactivated conda env fails to load MKL and takes
    the interpreter down with it.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    xc, yc = x - x.mean(), y - y.mean()
    denom = np.sqrt(np.square(xc).sum()) * np.sqrt(np.square(yc).sum())
    if denom == 0:
        return np.nan
    return float(np.multiply(xc, yc).sum() / denom)


def compare_metrics(outs_a: Path, outs_b: Path, label_a: str, label_b: str) -> pd.DataFrame:
    """Side-by-side metrics_summary.csv. CellRanger renamed some fields between
    9.0.1 and 10.1.0 (e.g. 'Valid UMIs' -> 'Valid UMI Sequences'), so the index is
    normalised before joining."""
    a, b = read_metrics(outs_a), read_metrics(outs_b)
    a.index = [_norm_metric(i) for i in a.index]
    b.index = [_norm_metric(i) for i in b.index]
    df = pd.DataFrame({label_a: a, label_b: b})
    numeric = df.apply(pd.to_numeric, errors="coerce")
    df["abs_diff"] = (numeric[label_b] - numeric[label_a]).abs()
    df["pct_diff"] = 100 * df["abs_diff"] / numeric[label_a].abs().replace(0, np.nan)
    return df


_METRIC_ALIASES = {
    "valid umis": "valid umi sequences",
}


def _norm_metric(name: str) -> str:
    key = str(name).strip().lower()
    return _METRIC_ALIASES.get(key, key)
