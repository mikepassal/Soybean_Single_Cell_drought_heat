"""Assembling the analysis object, from either source, with consistent obs columns.

Two sources:

* ``collab_h5ad``  -- ``adata_rna.h5ad``. Its ``uns['source']`` records that ``.X`` is the
  RNA ``data`` layer and ``layers['counts']`` is raw counts, so the counts are usable
  even though the object came out of their Seurat pipeline. This is the only source
  with all four conditions right now, so it is the default.
* ``realigned``    -- our own CellRanger output. Currently control-only; the QC/obs
  columns are built to match ``collab_h5ad`` exactly so scripts downstream do not care
  which one they were handed.

The realigned source keeps the organelle contigs (``GlmaCp*`` / ``GlmaCt*``), which the
collaborator dropped from the h5ad after computing ``pctCP``. That matters: chloroplast
load varies ~20x between libraries and tracks replicate batch, so it needs to be a
first-class covariate rather than something already averaged away.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from . import config as cfg
from . import matrices as mx

# Columns every downstream script relies on.
OBS_COLUMNS = [
    "treatment", "libraries", "replicate", "rep_block",
    "celltype_call", "cluster_annot", "pct_chloroplast",
]


def _rep_block(replicate: pd.Series) -> pd.Series:
    """rep1A is a re-run of the control rep1 well, so it folds into rep1.

    The point of this column is to have a batch factor that is *crossed* with
    treatment: rep1/rep2/rep3 each appear in all four conditions, whereas
    `libraries` is nested inside treatment. Anything nested inside treatment cannot
    be used to correct batch without also removing the treatment effect -- which is
    the core question these scripts exist to answer.
    """
    return pd.Categorical(
        replicate.astype(str).replace({"rep1A": "rep1"}),
        categories=["rep1", "rep2", "rep3"],
    )


def load_collab_h5ad(path: Path | None = None) -> ad.AnnData:
    path = path or cfg.COLLAB_H5AD
    adata = sc.read_h5ad(path)

    # .X is their log-normalised RNA data; we want raw counts as the working matrix
    # and will do our own normalisation, so the layer becomes .X.
    if "counts" not in adata.layers:
        raise ValueError(f"{path.name} has no layers['counts']")
    adata.X = adata.layers["counts"].copy()

    adata.obs["rep_block"] = _rep_block(adata.obs["replicate"])
    adata.obs["treatment"] = pd.Categorical(
        adata.obs["treatment"].astype(str), categories=cfg.TREATMENTS
    )
    # Organelle genes are gone from this object, but their summary survived.
    adata.obs["pct_chloroplast"] = adata.obs["pctCP"].astype(float)
    adata.obs["source"] = "collab_h5ad"
    adata.uns["counts_are_cellbender_corrected"] = True

    # The Seurat export leaves a uns['neighbors'] with no 'params', which scanpy's
    # Neighbors reads unconditionally and chokes on. Their graph is not something we
    # want to reuse anyway -- each variant builds its own.
    adata.uns.pop("neighbors", None)
    for key in list(adata.obsp):
        del adata.obsp[key]

    # Move their embeddings out of the default slots immediately. Otherwise the first
    # variant to run PCA overwrites obsm['X_pca'], and the 'theirs' variant silently
    # becomes a duplicate of 'unintegrated'.
    for key in ("X_pca", "X_umap"):
        if key in adata.obsm:
            adata.obsm[f"{key}_theirs"] = adata.obsm.pop(key)
    return adata


def load_realigned(samples: list[cfg.Sample] | None = None, min_genes: int = 500) -> ad.AnnData:
    """Concatenate our CellRanger matrices and attach the collaborator's labels.

    Cells that the collaborator dropped keep NaN cell-type labels rather than being
    discarded, so the effect of their QC is measurable instead of baked in.
    """
    samples = samples or cfg.available_samples("realigned")
    if not samples:
        raise RuntimeError("no realigned samples on disk yet")

    parts = []
    for s in samples:
        a = mx.read_matrix(s.matrix_dir("realigned"))
        a.obs["libraries"] = s.library
        a.obs["treatment"] = s.treatment
        a.obs["replicate"] = s.replicate
        a.obs["run_id"] = s.run_id
        parts.append(a)
        print(f"  {s.run_id:<12} {a.n_obs:>7,} cells x {a.n_vars:,} genes")

    adata = ad.concat(parts, join="outer", label=None, index_unique=None, merge="first")
    # Barcodes repeat across libraries; make them match the collaborator's naming.
    adata.obs_names = [
        cfg.barcode_to_collab_obs_name(bc.split("-")[0] + "-1", lib)
        for bc, lib in zip(adata.obs_names, adata.obs["libraries"])
    ]
    adata.obs_names_make_unique()

    adata = mx.annotate_organelles(adata)
    adata.obs["pct_chloroplast"] = adata.obs["pct_counts_chloroplast"].astype(float)
    adata.obs["pct_mito"] = adata.obs["pct_counts_mito"].astype(float)
    adata.obs["rep_block"] = _rep_block(adata.obs["replicate"])
    adata.obs["treatment"] = pd.Categorical(
        adata.obs["treatment"].astype(str), categories=cfg.TREATMENTS
    )
    adata.obs["source"] = "realigned"
    adata.uns["counts_are_cellbender_corrected"] = False

    _attach_collab_labels(adata)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    return adata


def _attach_collab_labels(adata: ad.AnnData) -> None:
    """Carry over celltype_call / cluster_annot for cells that survived their QC."""
    if not cfg.COLLAB_H5AD.exists():
        adata.obs["celltype_call"] = pd.NA
        adata.obs["cluster_annot"] = pd.NA
        return
    ref = sc.read_h5ad(cfg.COLLAB_H5AD, backed="r")
    for col in ("celltype_call", "cluster_annot"):
        s = pd.Series(ref.obs[col].astype(str).values, index=ref.obs_names)
        adata.obs[col] = pd.Categorical(s.reindex(adata.obs_names).values)
    kept = adata.obs["celltype_call"].notna().sum()
    print(f"  {kept:,} / {adata.n_obs:,} cells carry a collaborator cell-type label")


def load(source: str = "collab_h5ad", **kwargs) -> ad.AnnData:
    if source == "collab_h5ad":
        return load_collab_h5ad(**kwargs)
    if source == "realigned":
        return load_realigned(**kwargs)
    raise ValueError(f"unknown source: {source!r}")


# --------------------------------------------------------------------------- QC


def qc(
    adata: ad.AnnData,
    min_genes: int = 300,
    max_pct_chloroplast: float | None = None,
    drop_organelle_genes: bool = True,
) -> ad.AnnData:
    """Light QC. Deliberately lighter than the collaborator's, so that their filtering
    is a variable we can vary rather than a fixed part of the input.

    ``drop_organelle_genes`` removes chloroplast/mito genes *from the matrix* after
    per-cell fractions have been recorded. This is the standard move and it is what
    the collaborator did -- but note it does not undo the damage: the organelle reads
    already consumed sequencing depth, so a high-chloroplast cell still has fewer
    nuclear UMIs. That residual is what ``pct_chloroplast`` is for downstream.
    """
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)

    if max_pct_chloroplast is not None and "pct_chloroplast" in adata.obs:
        before = adata.n_obs
        adata = adata[adata.obs["pct_chloroplast"] <= max_pct_chloroplast].copy()
        print(f"  chloroplast filter <={max_pct_chloroplast}%: {before:,} -> {adata.n_obs:,} cells")

    if drop_organelle_genes and "chloroplast" in adata.var:
        keep = ~(adata.var["chloroplast"].to_numpy() | adata.var["mito"].to_numpy())
        print(f"  dropping {(~keep).sum()} organelle genes")
        adata = adata[:, keep].copy()

    sc.pp.filter_genes(adata, min_cells=5)
    adata.layers["counts"] = adata.X.copy()
    return adata


def normalize(adata: ad.AnnData, n_top_genes: int = 2000, batch_key: str | None = "libraries") -> ad.AnnData:
    """Standard log-normalisation + HVG selection, keeping raw counts in a layer."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, batch_key=batch_key if batch_key in adata.obs else None
    )
    return adata
