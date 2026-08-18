"""Q1: do we reproduce the collaborator's count matrix?

Compares, per library and at three stages:

1. ``metrics_summary.csv``  -- the CellRanger run itself (mapping, cell calling)
2. ``filtered_feature_bc_matrix`` / ``raw_feature_bc_matrix`` -- our realignment vs theirs
3. their CellBender output and ``adata_rna.h5ad`` -- how much the *post*-CellRanger
   steps move the counts, which sets the scale for judging any difference in (2)

Stage 3 matters because their analysis object is not CellRanger output: it is
CellRanger -> CellBender -> Seurat QC. A perfect match at stage 2 still leaves
CellBender as a source of divergence, so the report quantifies both.

Runs on whatever is on disk. Right now that is the four control libraries.

    python Reprocess_and_recluster/scripts/01_compare_count_matrices.py
    python Reprocess_and_recluster/scripts/01_compare_count_matrices.py --conditions control heat
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as cfg
from pipeline import matrices as mx


def compare_sample(sample: cfg.Sample, do_cellbender: bool) -> dict:
    rows: dict[str, dict] = {}

    print(f"\n=== {sample.run_id}  ({sample.library}) ===", flush=True)

    # --- stage 1: CellRanger run metrics ----------------------------------
    metrics = None
    if (sample.realigned_outs / "metrics_summary.csv").exists() and (
        sample.collab_outs / "metrics_summary.csv"
    ).exists():
        metrics = mx.compare_metrics(
            sample.collab_outs, sample.realigned_outs, "collaborator", "realigned"
        )
        n_diff = int((metrics["abs_diff"].fillna(0) > 0).sum())
        print(f"  metrics_summary: {n_diff} of {len(metrics)} fields differ")

    # --- stage 2: the matrices themselves ---------------------------------
    for filtered in (True, False):
        stage = "filtered" if filtered else "raw"
        mine_dir = sample.matrix_dir("realigned", filtered)
        theirs_dir = sample.matrix_dir("collaborator", filtered)
        if not (mine_dir / "matrix.mtx.gz").exists() or not (theirs_dir / "matrix.mtx.gz").exists():
            print(f"  {stage}: missing on one side, skipped")
            continue
        mine = mx.read_matrix(mine_dir)
        theirs = mx.read_matrix(theirs_dir)
        res = mx.compare_matrices(theirs, mine, "collaborator", "realigned")
        rows[f"cellranger_{stage}"] = res
        print(
            f"  {stage:8s}: cells {res['n_cells_collaborator']} vs {res['n_cells_realigned']}"
            f" | shared {res['n_cells_shared']}"
            f" | identical={res.get('identical')}"
            f" | differing entries={res.get('n_differing_entries', 0):,}"
        )
        if filtered:
            mine_filtered, theirs_filtered = mine, theirs

    # --- stage 3: what CellBender and their QC did on top ------------------
    if do_cellbender and sample.cellbender_h5.exists():
        cb = mx.read_cellbender(sample.cellbender_h5)
        # CellBender writes gene ids in var_names for 10x-style h5; align on ids.
        res = mx.compare_matrices(theirs_filtered, cb, "collaborator", "cellbender")
        rows["cellranger_vs_cellbender"] = res
        print(
            f"  cellbender: calls {res['n_cells_cellbender']:,} cells vs CellRanger's"
            f" {res['n_cells_collaborator']:,} ({res['n_cells_shared']:,} shared);"
            f" on the shared cells it removes"
            f" {100 * res['frac_umis_removed_shared_cells']:.1f}% of UMIs"
        )
    elif do_cellbender:
        print(f"  cellbender: not found at {sample.cellbender_h5}")

    return {"metrics": metrics, "matrix": rows}


def realigned_metrics_table(samples: list[cfg.Sample]) -> pd.DataFrame:
    """Our own CellRanger metrics for every library, whether or not theirs is on disk.

    The collaborator's ``metrics_summary.csv`` is not part of every copy of their
    data, but read depth is the thing that explains the rep2 libraries, so record
    ours unconditionally.
    """
    rows = {}
    for s in samples:
        path = s.realigned_outs / "metrics_summary.csv"
        if path.exists():
            rows[s.run_id] = mx.read_metrics(s.realigned_outs)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).T
    df.index.name = "run_id"
    df.insert(0, "replicate", [cfg.SAMPLES_BY_RUN_ID[r].replicate for r in df.index])
    return df


def compare_to_analysis_object(samples: list[cfg.Sample]) -> pd.DataFrame:
    """How our realigned counts line up with the cells that survived into adata_rna.h5ad."""
    if not cfg.COLLAB_H5AD.exists():
        return pd.DataFrame()
    print(f"\n=== vs {cfg.COLLAB_H5AD.name} (their final analysis object) ===", flush=True)
    final = ad.read_h5ad(cfg.COLLAB_H5AD, backed="r")
    rows = []
    for s in samples:
        keep = final.obs_names[final.obs["libraries"].astype(str) == s.library]
        their_bcs = pd.Index([b.rsplit("_", 1)[-1] for b in keep])
        mine = mx.read_matrix(sample_dir := s.matrix_dir("realigned"))
        shared = mine.obs_names.intersection(their_bcs)
        rows.append(
            {
                "run_id": s.run_id,
                "library": s.library,
                "cells_realigned_filtered": mine.n_obs,
                "cells_in_adata_rna": len(their_bcs),
                "shared_barcodes": len(shared),
                "frac_of_adata_rna_recovered": len(shared) / max(len(their_bcs), 1),
                "frac_of_ours_retained": len(shared) / max(mine.n_obs, 1),
            }
        )
        print(
            f"  {s.run_id:<12} ours {mine.n_obs:>6,} | theirs {len(their_bcs):>6,}"
            f" | shared {len(shared):>6,}"
            f" ({100 * len(shared) / max(len(their_bcs), 1):.1f}% of theirs)"
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conditions", nargs="*", default=None, help="restrict to these conditions")
    ap.add_argument("--no-cellbender", action="store_true", help="skip the CellBender comparison")
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "count_matrix_comparison")
    args = ap.parse_args()

    print(cfg.describe_availability("realigned"))

    samples = cfg.available_samples("realigned")
    if args.conditions:
        samples = [s for s in samples if s.condition in args.conditions]
    samples = [s for s in samples if s.has("collaborator")]
    if not samples:
        print("\nNothing to compare: no library has both a realigned and a collaborator matrix.")
        return

    args.outdir.mkdir(parents=True, exist_ok=True)

    metric_frames, matrix_rows = {}, []
    for s in samples:
        res = compare_sample(s, do_cellbender=not args.no_cellbender)
        if res["metrics"] is not None:
            metric_frames[s.run_id] = res["metrics"]
        for stage, row in res["matrix"].items():
            matrix_rows.append({"run_id": s.run_id, "library": s.library, "stage": stage, **row})

    if metric_frames:
        allm = pd.concat(metric_frames, names=["run_id", "metric"])
        allm.to_csv(args.outdir / "metrics_summary_comparison.csv")
        print(f"\nwrote {args.outdir / 'metrics_summary_comparison.csv'}")
        disagree = allm[allm["abs_diff"].fillna(0) > 0]
        print(f"metrics fields differing anywhere: {len(disagree)}")
        if len(disagree):
            print(disagree.to_string())

    if matrix_rows:
        mdf = pd.DataFrame(matrix_rows)
        mdf.to_csv(args.outdir / "matrix_comparison.csv", index=False)
        print(f"wrote {args.outdir / 'matrix_comparison.csv'}")
        cols = [c for c in ["run_id", "stage", "n_cells_shared", "identical",
                            "n_differing_entries", "frac_entries_differing",
                            "per_gene_umi_pearson"] if c in mdf]
        print("\n" + mdf[cols].to_string(index=False))

    rmetrics = realigned_metrics_table(samples)
    if len(rmetrics):
        rmetrics.to_csv(args.outdir / "realigned_metrics_summary.csv")
        print(f"wrote {args.outdir / 'realigned_metrics_summary.csv'}")
        depth = rmetrics[["replicate", "Number of Reads", "Estimated Number of Cells",
                          "Sequencing Saturation"]]
        print("\n" + depth.to_string())

    final_df = compare_to_analysis_object(samples)
    if len(final_df):
        final_df.to_csv(args.outdir / "vs_adata_rna.csv", index=False)
        print(f"\nwrote {args.outdir / 'vs_adata_rna.csv'}")

    missing = cfg.missing_samples("realigned")
    if missing:
        print(f"\nStill awaiting realignment ({len(missing)}): "
              + ", ".join(s.run_id for s in missing))
        print("Rerun this script unchanged once they appear.")


if __name__ == "__main__":
    main()
