"""Differential expression, crossed over grouping and over what the test consumes.

The grid has two axes.

**What the test consumes**

* ``pseudobulk_counts`` -- sum *raw counts* per library, DESeq2 negative binomial.
  The reference analysis. Replicates are libraries, i.e. real biological replicates.
  Note this is *invariant to Harmony*: Harmony only moves an embedding, so it can
  reach DE solely by changing which cells are grouped together. If Harmony-clustered
  and unintegrated-clustered results differ, the cause is cluster assignment, not
  correction of the counts.
* ``pseudobulk_corrected`` -- average *ComBat-corrected expression* per library, then
  a linear model. This is the analogue of running DE on a Seurat integrated assay,
  and it is where over-correction actually destroys signal.
* ``wilcoxon_cells`` -- cell-level rank-sum on log-normalised values, treating cells
  as replicates. What ``FindMarkers`` does by default. Included because it is the
  usual comparator, and because its p-values are severely anti-conservative here
  (n = thousands of cells, but only 3-4 independent plants per condition).

**How cells are grouped**

* ``library``            -- one pseudosample per library; no cell-type split.
* ``celltype_call``      -- the collaborator's labels.
* ``cluster_annot``      -- the collaborator's finer labels.
* ``leiden_<variant>``   -- de-novo clusters from each integration variant.

Covariate handling: ``~ rep_block + treatment`` is offered alongside ``~ treatment``
because chloroplast load, and hence library composition, tracks replicate batch far
more strongly than it tracks condition.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from . import config as cfg

ALPHA, LFC_CUT = 0.05, 1.0
MIN_CELLS_PER_PSEUDOSAMPLE = 10
MIN_REPS_PER_TREATMENT = 2

# The effect-size column means different things per arm, so the cut has to vary with it.
# DESeq2 arms report a genuine log2 fold change; the corrected arm reports a difference in
# mean log1p(CP10K), which is an order of magnitude smaller and not comparable.
LFC_CUT_BY_TEST = {
    "pseudobulk_counts": 1.0,
    "pseudobulk_counts_repblock": 1.0,
    "pseudobulk_corrected": 0.0,   # significance rests on padj alone
    "wilcoxon_cells": 1.0,
}


def lfc_cut_for(test: str) -> float:
    return LFC_CUT_BY_TEST.get(test, LFC_CUT)


# --------------------------------------------------------------------------- pseudobulk


def pseudobulk(
    adata: ad.AnnData,
    group_keys: list[str],
    layer: str | None = "counts",
    how: str = "sum",
    meta_keys: tuple[str, ...] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate cells into pseudosamples.

    ``how='sum'`` for raw counts (what DESeq2 wants); ``how='mean'`` for corrected
    expression values, where summing would just re-encode cell number.
    """
    labels = adata.obs[group_keys].astype(str).agg(" | ".join, axis=1).to_numpy()
    cats = pd.Categorical(labels)
    n_cells = np.bincount(cats.codes, minlength=len(cats.categories))

    X = adata.layers[layer] if layer else adata.X
    indicator = sparse.csr_matrix(
        (np.ones(adata.n_obs), (cats.codes, np.arange(adata.n_obs))),
        shape=(len(cats.categories), adata.n_obs),
    )
    agg = indicator @ X
    agg = np.asarray(agg.todense()) if sparse.issparse(agg) else np.asarray(agg)
    if how == "mean":
        agg = agg / np.maximum(n_cells, 1)[:, None]
    else:
        agg = np.rint(agg).astype(np.int64)

    counts = pd.DataFrame(agg, index=list(cats.categories), columns=adata.var_names)

    meta_cols = list(dict.fromkeys(group_keys + [k for k in meta_keys if k in adata.obs]))
    frame = adata.obs[meta_cols].astype(str).copy()
    frame["_group"] = labels
    meta = frame.groupby("_group", observed=True)[meta_cols].agg(lambda s: s.value_counts().index[0])
    meta["n_cells"] = pd.Series(n_cells, index=list(cats.categories))
    meta = meta.loc[counts.index]
    if "treatment" in meta:
        meta["treatment"] = pd.Categorical(meta["treatment"], categories=cfg.TREATMENTS)
    return counts, meta


def filter_genes(counts: pd.DataFrame, min_total: int = 10, min_frac_samples: float = 0.5) -> pd.DataFrame:
    min_samples = max(2, int(np.ceil(min_frac_samples * counts.shape[0])))
    keep = ((counts > 0).sum(axis=0) >= min_samples) & (counts.sum(axis=0) >= min_total)
    return counts.loc[:, keep.to_numpy()]


# --------------------------------------------------------------------------- tests


def deseq(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    design: str = "~treatment",
    contrasts: dict | None = None,
    n_cpus: int = 8,
) -> dict[str, pd.DataFrame]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    contrasts = contrasts or cfg.CONTRASTS
    meta = meta.copy()
    for col in ("treatment", "rep_block"):
        if col in meta:
            meta[col] = meta[col].astype(str)

    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, n_cpus=n_cpus, quiet=True)
    dds.deseq2()
    out = {}
    for name, contrast in contrasts.items():
        if contrast[1] not in set(meta["treatment"]):
            continue
        stat = DeseqStats(dds, contrast=contrast, n_cpus=n_cpus, quiet=True)
        stat.summary()
        out[name] = stat.results_df.sort_values("padj")
    return out


def linear_model(
    values: pd.DataFrame,
    meta: pd.DataFrame,
    covariates: list[str] | None = None,
    contrasts: dict | None = None,
) -> dict[str, pd.DataFrame]:
    """OLS per gene on already-normalised (e.g. ComBat-corrected) pseudobulk values.

    DESeq2's negative-binomial model is wrong for corrected expression -- the values
    are continuous and can be negative -- so this is the right test for the
    ``pseudobulk_corrected`` arm, and it is what makes that arm comparable to a Seurat
    integrated-assay analysis.

    IMPORTANT: the ``log2FoldChange`` column is named for interface compatibility with
    the DESeq2 results tables, but it is **not** a log2 fold change. It is a difference
    in mean ``log1p(CP10K)`` between groups, which lives on a much smaller scale --
    typically |coef| < 0.35 at the 99th percentile. Applying the DESeq2-appropriate
    ``|LFC| >= 1`` cut here would reject essentially everything for the wrong reason,
    so ``LFC_CUT_BY_TEST`` sets this arm's cut to 0 and significance rests on padj.
    """
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    contrasts = contrasts or cfg.CONTRASTS
    meta = meta.copy()
    treat = pd.Categorical(meta["treatment"].astype(str), categories=cfg.TREATMENTS)
    present = [t for t in cfg.TREATMENTS if t in set(treat.dropna())]
    design = pd.get_dummies(
        pd.Categorical(treat, categories=present), prefix="treatment", drop_first=True
    ).astype(float)
    design.index = meta.index
    for cov in covariates or []:
        if cov in meta:
            d = pd.get_dummies(meta[cov].astype(str), prefix=cov, drop_first=True).astype(float)
            design = pd.concat([design, d.set_index(design.index)], axis=1)
    X = sm.add_constant(design, has_constant="add")

    Y = values.loc[meta.index]
    results: dict[str, list] = {name: [] for name in contrasts}
    coefs, ses, dfr = _ols_batch(X.to_numpy(dtype=float), Y.to_numpy(dtype=float))
    colidx = {c: i for i, c in enumerate(X.columns)}

    from scipy import stats as sps

    for name, (_, level, _ref) in contrasts.items():
        col = f"treatment_{level}"
        if col not in colidx:
            continue
        j = colidx[col]
        beta, se = coefs[:, j], ses[:, j]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = beta / se
        p = 2 * sps.t.sf(np.abs(t), dfr)
        ok = np.isfinite(p)
        padj = np.full_like(p, np.nan)
        if ok.any():
            padj[ok] = multipletests(p[ok], method="fdr_bh")[1]
        results[name] = pd.DataFrame(
            {"baseMean": Y.mean(axis=0).to_numpy(), "log2FoldChange": beta,
             "lfcSE": se, "stat": t, "pvalue": p, "padj": padj},
            index=Y.columns,
        ).sort_values("padj")
    return {k: v for k, v in results.items() if len(v)}


def _ols_batch(X: np.ndarray, Y: np.ndarray):
    """Per-gene OLS in one solve. Y is samples x genes."""
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ Y            # p x genes
    resid = Y - X @ beta
    dfr = max(n - np.linalg.matrix_rank(X), 1)
    sigma2 = (resid ** 2).sum(axis=0) / dfr
    se = np.sqrt(np.outer(np.diag(XtX_inv), sigma2))  # p x genes
    return beta.T, se.T, dfr


def wilcoxon_cells(
    adata: ad.AnnData, layer: str = "lognorm", contrasts: dict | None = None
) -> dict[str, pd.DataFrame]:
    """Cell-level rank-sum, Seurat FindMarkers style.

    Included as a comparator, not as a recommendation: cells within a library are not
    independent, so these p-values are wildly anti-conservative.
    """
    import scanpy as sc

    contrasts = contrasts or cfg.CONTRASTS
    sub = adata.copy()
    sub.X = sub.layers[layer]
    out = {}
    for name, (_, level, ref) in contrasts.items():
        mask = sub.obs["treatment"].astype(str).isin([level, ref])
        if mask.sum() == 0 or sub.obs.loc[mask, "treatment"].nunique() < 2:
            continue
        s = sub[mask].copy()
        s.obs["_g"] = pd.Categorical(s.obs["treatment"].astype(str), categories=[ref, level])
        sc.tl.rank_genes_groups(s, "_g", groups=[level], reference=ref, method="wilcoxon")
        df = sc.get.rank_genes_groups_df(s, group=level)
        out[name] = pd.DataFrame(
            {"log2FoldChange": df["logfoldchanges"].to_numpy(),
             "stat": df["scores"].to_numpy(),
             "pvalue": df["pvals"].to_numpy(),
             "padj": df["pvals_adj"].to_numpy()},
            index=df["names"],
        ).sort_values("padj")
    return out


# --------------------------------------------------------------------------- summary


def sig_genes(res: pd.DataFrame, alpha: float = ALPHA, lfc_cut: float = LFC_CUT) -> pd.Index:
    return res.index[(res["padj"] < alpha) & (res["log2FoldChange"].abs() >= lfc_cut)]


def summarize(results: dict[str, pd.DataFrame], **tags) -> pd.DataFrame:
    lfc_cut = lfc_cut_for(str(tags.get("test", "")))
    rows = []
    for name, res in results.items():
        sig = res.loc[sig_genes(res, lfc_cut=lfc_cut)]
        rows.append({
            **tags,
            "lfc_cut": lfc_cut,
            "comparison": name,
            "genes_tested": int(res["padj"].notna().sum()),
            "n_padj_sig": int((res["padj"] < ALPHA).sum()),
            "n_sig": len(sig),
            "n_up": int((sig["log2FoldChange"] > 0).sum()),
            "n_down": int((sig["log2FoldChange"] < 0).sum()),
        })
    return pd.DataFrame(rows)


def save(results: dict[str, pd.DataFrame], outdir: Path, prefix: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for name, res in results.items():
        res.to_csv(outdir / f"{prefix}__{name.replace(' ', '_')}.csv")
