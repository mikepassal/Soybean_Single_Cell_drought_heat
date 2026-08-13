"""How batch-confounded is each clustering, and what does that cost the DE?

Motivated by a concrete observation in the collaborator's `cluster_annot`: several
clusters (Mesophyll-7, Mesophyll-8, Mesophyll-13, Epidermal/Pavement-9) are drawn almost
entirely from the rep3 libraries. Mesophyll-7, for instance, is 725 control_rep3 + 473
drought_rep3 + 632 Heat_rep3 + 349 HD_rep3, against 1-22 cells from every rep1/rep2
library.

Clusters like that are replicate-batch artifacts, not cell types, and they are expensive
in a per-cell-type DE: each condition contributes only its single rep3 library, so the
cluster has one replicate per treatment and the contrast is unestimable. The stress
conditions silently drop out of the analysis rather than returning a null result.

For every grouping scheme this reports:

* ``max_repblock_frac``   -- the largest share any one replicate block holds. ~1.0 means
  the cluster is a batch, not a cell type.
* ``effective_libraries`` -- libraries contributing at least ``--min-cells`` cells.
* ``treatments_testable`` -- treatments with at least 2 such libraries. Below 2 for
  Control, or below 1 for the stresses, and the cluster contributes nothing.
* ``chloroplast_median``  -- since chloroplast load is what most strongly separates the
  replicate blocks.

    python Reprocess_and_recluster/scripts/05_cluster_composition.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as cfg


def composition(adata, group_col: str, min_cells: int) -> pd.DataFrame:
    obs = adata.obs
    rows = []
    for ct, sub in obs.groupby(group_col, observed=True):
        by_lib = sub.groupby("libraries", observed=True).size()
        eff = by_lib[by_lib >= min_cells]
        lib_treat = obs.groupby("libraries", observed=True)["treatment"].agg(
            lambda s: s.value_counts().index[0]
        )
        treat_counts = pd.Series(
            [lib_treat.get(l) for l in eff.index], dtype=object
        ).value_counts()
        testable = [t for t in cfg.TREATMENTS if treat_counts.get(t, 0) >= 2]

        rb = sub["rep_block"].value_counts(normalize=True)
        rows.append({
            "grouping": group_col,
            "cluster": str(ct),
            "n_cells": len(sub),
            "max_repblock_frac": float(rb.max()) if len(rb) else np.nan,
            "dominant_repblock": str(rb.idxmax()) if len(rb) else "-",
            "effective_libraries": int(len(eff)),
            "treatments_testable": len(testable),
            "testable": ",".join(testable) if testable else "-",
            "chloroplast_median": float(sub["pct_chloroplast"].median()),
        })
    return pd.DataFrame(rows).sort_values("max_repblock_frac", ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path,
                    default=cfg.RESULTS_ROOT / "integration" / "integrated_collab_h5ad.h5ad")
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "integration")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    groupings = ["celltype_call", "cluster_annot"] + [
        c for c in adata.obs.columns if c.startswith("leiden_")
    ]

    frames = [composition(adata, g, args.min_cells) for g in groupings if g in adata.obs]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(args.outdir / "cluster_composition.csv", index=False)

    print("=== per-grouping summary ===")
    summary = df.groupby("grouping").agg(
        n_clusters=("cluster", "size"),
        median_max_repblock=("max_repblock_frac", "median"),
        n_batch_dominated=("max_repblock_frac", lambda s: int((s > 0.8).sum())),
        n_all4_testable=("treatments_testable", lambda s: int((s == 4).sum())),
        n_unusable=("treatments_testable", lambda s: int((s < 2).sum())),
        cells_in_batch_dominated=("n_cells", lambda s: 0),
    )
    # cells sitting inside a batch-dominated cluster, per grouping
    dom = df[df["max_repblock_frac"] > 0.8].groupby("grouping")["n_cells"].sum()
    total = df.groupby("grouping")["n_cells"].sum()
    summary["pct_cells_batch_dominated"] = (100 * dom / total).reindex(summary.index).fillna(0).round(1)
    summary = summary.drop(columns=["cells_in_batch_dominated"])
    print(summary.to_string())

    print("\n=== the worst offenders (max_repblock_frac > 0.8) ===")
    worst = df[df["max_repblock_frac"] > 0.8].sort_values(
        ["grouping", "max_repblock_frac"], ascending=[True, False]
    )
    if len(worst):
        print(worst[["grouping", "cluster", "n_cells", "max_repblock_frac", "dominant_repblock",
                     "effective_libraries", "treatments_testable", "chloroplast_median"]]
              .to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    else:
        print("none")

    print(f"\nwrote {args.outdir / 'cluster_composition.csv'}")
    print("\nRead: max_repblock_frac ~1.0 means the cluster is a replicate batch, not a")
    print("cell type. treatments_testable < 2 means it contributes nothing to per-cell-type DE.")


if __name__ == "__main__":
    main()
