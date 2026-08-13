"""Integration variants, chosen to isolate *why* integration could cost DE signal.

The variants differ along two axes that matter independently:

**What the batch key is.** ``libraries`` is nested inside ``treatment`` -- every cell in
a library shares a treatment. Correcting on it asks the method to remove variation that
is perfectly confounded with the biology, so it removes treatment signal by
construction. ``rep_block`` (rep1/rep2/rep3, with rep1A folded into rep1) is *crossed*
with treatment: each block contains all four conditions, so correcting on it removes
batch while leaving the treatment contrast estimable. This is the axis most likely to
explain a weak Heat/Drought result.

**What the method corrects.** Harmony adjusts an *embedding* only -- the expression
matrix is untouched, so raw-count pseudobulk DE is mathematically unaffected by it and
can only change through cluster assignment. ComBat rewrites the *expression values*, so
it can and does propagate into any DE run on corrected values. Keeping both makes the
distinction measurable rather than assumed.

The collaborator's own result is carried as the ``theirs`` variant: their PCA and their
``cluster_annot`` labels, read straight out of the h5ad.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


@dataclass(frozen=True)
class Variant:
    name: str
    method: str  # none | harmony | combat | theirs
    batch_key: str | None = None
    notes: str = ""
    # Set when the method rewrites expression values rather than just an embedding.
    corrects_expression: bool = False

    @property
    def rep_key(self) -> str:
        return "X_pca" if self.method in ("none",) else f"X_{self.name}"

    @property
    def cluster_key(self) -> str:
        return f"leiden_{self.name}"


VARIANTS: list[Variant] = [
    Variant(
        "unintegrated", "none", None,
        "PCA on log-normalised HVGs. Baseline: no correction of any kind.",
    ),
    Variant(
        "harmony_repblock", "harmony", "rep_block",
        "Harmony on a batch factor crossed with treatment. The defensible choice: "
        "removes replicate batch while leaving the treatment contrast estimable.",
    ),
    Variant(
        "harmony_library", "harmony", "libraries",
        "Harmony on a batch factor nested inside treatment. Expected to suppress "
        "condition separation, since library and treatment are confounded.",
    ),
    Variant(
        "combat_repblock", "combat", "rep_block",
        "ComBat on the crossed batch factor; rewrites expression values.",
        corrects_expression=True,
    ),
    Variant(
        "combat_library", "combat", "libraries",
        "ComBat on the nested batch factor; rewrites expression values. This is the "
        "closest analogue to running DE on a Seurat integrated assay.",
        corrects_expression=True,
    ),
    Variant(
        "theirs", "theirs", None,
        "The collaborator's own PCA and cluster_annot labels, from adata_rna.h5ad.",
    ),
]

VARIANTS_BY_NAME = {v.name: v for v in VARIANTS}


def run_variant(
    adata: ad.AnnData,
    variant: Variant,
    n_pcs: int = 30,
    resolution: float = 0.8,
    seed: int = 0,
) -> ad.AnnData:
    """Compute the variant's embedding and its de-novo leiden clustering, in place.

    Adds ``obsm[variant.rep_key]`` and ``obs[variant.cluster_key]``; for
    expression-correcting methods also adds ``layers[f'{name}_corrected']``.
    """
    print(f"\n--- {variant.name} ({variant.method}"
          + (f", batch={variant.batch_key}" if variant.batch_key else "") + ") ---", flush=True)

    if variant.method == "theirs":
        # Must be their PCA, not whatever is currently in obsm['X_pca'] -- build.load
        # parks theirs under X_pca_theirs precisely so this cannot pick up ours.
        if "X_pca_theirs" not in adata.obsm:
            raise ValueError("adata has no X_pca_theirs to reuse for the 'theirs' variant")
        their_pca = np.asarray(adata.obsm["X_pca_theirs"])
        adata.obsm[variant.rep_key] = their_pca[:, : min(n_pcs, their_pca.shape[1])]

    elif variant.method == "none":
        _pca(adata, n_pcs, seed)
        adata.obsm[variant.rep_key] = adata.obsm["X_pca"][:, :n_pcs]

    elif variant.method == "harmony":
        _pca(adata, n_pcs, seed)
        import harmonypy

        ho = harmonypy.run_harmony(
            adata.obsm["X_pca"][:, :n_pcs],
            adata.obs,
            vars_use=[variant.batch_key],
            max_iter_harmony=20,
        )
        # harmonypy <2 returns Z_corr as PCs x cells, >=2.0 as cells x PCs.
        Z = np.asarray(ho.Z_corr)
        adata.obsm[variant.rep_key] = Z if Z.shape[0] == adata.n_obs else Z.T

    elif variant.method == "combat":
        # ComBat needs dense log-normalised values; restrict to HVGs to keep it
        # tractable, then PCA the corrected block.
        hvg = adata.var["highly_variable"].to_numpy()
        sub = adata[:, hvg].copy()
        sub.X = sub.layers["lognorm"].copy()
        sc.pp.combat(sub, key=variant.batch_key)
        layer = f"{variant.name}_corrected"
        adata.layers[layer] = _scatter_back(adata, sub, hvg)
        sc.pp.scale(sub, max_value=10)
        sc.tl.pca(sub, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
        adata.obsm[variant.rep_key] = sub.obsm["X_pca"]

    else:
        raise ValueError(f"unknown method {variant.method!r}")

    sc.pp.neighbors(adata, use_rep=variant.rep_key, key_added=variant.name, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=variant.cluster_key,
        neighbors_key=variant.name,
        flavor="igraph",
        n_iterations=2,
        random_state=seed,
    )
    n_clusters = adata.obs[variant.cluster_key].nunique()
    print(f"  {n_clusters} de-novo clusters")
    return adata


def _pca(adata: ad.AnnData, n_pcs: int, seed: int) -> None:
    if "X_pca_ours" in adata.obsm:
        adata.obsm["X_pca"] = adata.obsm["X_pca_ours"]
        return
    hvg = adata.var["highly_variable"].to_numpy()
    sub = adata[:, hvg].copy()
    sub.X = sub.layers["lognorm"].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = sub.obsm["X_pca"]
    adata.obsm["X_pca_ours"] = sub.obsm["X_pca"]


def _scatter_back(adata: ad.AnnData, sub: ad.AnnData, hvg_mask: np.ndarray) -> np.ndarray:
    """Put corrected HVG values back into a full-width matrix, leaving non-HVGs as
    their uncorrected log-normalised values."""
    full = np.asarray(
        adata.layers["lognorm"].todense() if hasattr(adata.layers["lognorm"], "todense")
        else adata.layers["lognorm"]
    ).copy()
    full[:, hvg_mask] = np.asarray(sub.X)
    return full


# --------------------------------------------------------------------------- diagnostics


def batch_mixing(adata: ad.AnnData, variant: Variant, key: str, n_neighbors: int = 30) -> float:
    """Mean fraction of a cell's neighbours sharing its label under ``key``.

    Run with ``key='rep_block'`` this measures batch removal (lower is better mixing);
    with ``key='treatment'`` it measures how much condition separation survived
    (higher means the biology is still there). Reporting both is the point -- a method
    that flattens both is over-correcting.
    """
    from sklearn.neighbors import NearestNeighbors

    rep = adata.obsm[variant.rep_key]
    labels = adata.obs[key].astype(str).to_numpy()
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(rep)
    _, idx = nn.kneighbors(rep)
    same = (labels[idx[:, 1:]] == labels[:, None]).mean(axis=1)
    return float(same.mean())


def variant_diagnostics(adata: ad.AnnData, variants: list[Variant]) -> pd.DataFrame:
    rows = []
    for v in variants:
        if v.rep_key not in adata.obsm:
            continue
        row = {
            "variant": v.name,
            "method": v.method,
            "batch_key": v.batch_key or "-",
            "n_clusters": int(adata.obs[v.cluster_key].nunique()) if v.cluster_key in adata.obs else np.nan,
            "batch_purity_repblock": batch_mixing(adata, v, "rep_block"),
            "treatment_purity": batch_mixing(adata, v, "treatment"),
        }
        if "celltype_call" in adata.obs and v.cluster_key in adata.obs:
            row["ARI_vs_celltype_call"] = _ari(adata.obs["celltype_call"], adata.obs[v.cluster_key])
        rows.append(row)
    df = pd.DataFrame(rows)
    # Null expectation for purity: the share each label would get by chance.
    for key, col in [("rep_block", "batch_purity_repblock"), ("treatment", "treatment_purity")]:
        p = adata.obs[key].value_counts(normalize=True)
        df[col + "_null"] = float((p ** 2).sum())
    return df


def _ari(a: pd.Series, b: pd.Series) -> float:
    from sklearn.metrics import adjusted_rand_score

    mask = a.notna() & b.notna()
    return float(adjusted_rand_score(a[mask].astype(str), b[mask].astype(str)))
