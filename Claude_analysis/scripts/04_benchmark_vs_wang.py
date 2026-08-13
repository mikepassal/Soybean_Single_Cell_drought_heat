"""Score every DE result against the Wang et al. 24h heat / drought / combined lists.

This is the arbiter. Comparing strategies by "how many genes did it call" rewards
whichever one is most anti-conservative -- cell-level Wilcoxon will always win, and it
is the least trustworthy. Scoring against an independent bulk RNA-seq measurement of
the same three stresses asks a better question: which strategy recovers real biology,
and at what enrichment.

Reads the tables written by 03_run_de.py and writes one row per result file.

    python Reprocess_and_recluster/scripts/04_benchmark_vs_wang.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as cfg
from pipeline import de as de_mod
from pipeline import reference as ref


def parse_name(path: Path) -> dict:
    """`pseudobulk_counts__celltype_call__Mesophyll__Heat_vs_Control.csv` -> parts."""
    stem = path.stem
    if "__" not in stem:
        return {}
    *head, contrast = stem.split("__")
    contrast = contrast.replace("_", " ")
    test = head[0] if head else "?"
    grouping = head[1] if len(head) > 1 else "library"
    cell_type = head[2] if len(head) > 2 else "all"
    return {"test": test, "grouping": grouping, "cell_type": cell_type, "comparison": contrast}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", type=Path, default=cfg.RESULTS_ROOT / "de" / "tables")
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "de")
    ap.add_argument("--min-probability", type=float, default=0.0,
                    help="restrict Wang genes to at least this posterior probability")
    ap.add_argument("--alpha", type=float, default=de_mod.ALPHA)
    ap.add_argument("--lfc", type=float, default=None,
                    help="override the per-test effect-size cut")
    args = ap.parse_args()

    files = sorted(args.tables.glob("*.csv"))
    if not files:
        raise SystemExit(f"no result tables in {args.tables}\nRun 03_run_de.py first.")
    print(f"scoring {len(files)} result tables against Wang et al.\n")

    rows = []
    for path in files:
        meta = parse_name(path)
        if meta.get("comparison") not in ref.WANG_SETS:
            continue
        res = pd.read_csv(path, index_col=0)
        if "padj" not in res:
            continue
        tested = res.index[res["padj"].notna()]
        # The corrected arm's effect size is not a log2 fold change, so it gets its own
        # cut -- scoring it at |LFC| >= 1 would call nothing regardless of the biology.
        lfc = args.lfc if args.lfc is not None else de_mod.lfc_cut_for(meta.get("test", ""))
        called = de_mod.sig_genes(res, args.alpha, lfc)
        stats = ref.enrichment(called, tested, meta["comparison"], args.min_probability)
        rows.append({**meta, "n_tested": len(tested), **stats, "file": path.name})

    if not rows:
        raise SystemExit("nothing scored")

    df = pd.DataFrame(rows).drop(columns=["contrast"])
    df.to_csv(args.outdir / "wang_benchmark.csv", index=False)
    print(f"wrote {args.outdir / 'wang_benchmark.csv'}  ({len(df)} rows)\n")

    # Library-level: the honest head-to-head, no cell-type splitting.
    lib = df[df["grouping"] == "library"]
    if len(lib):
        print("=== library-level pseudobulk, scored against Wang ===")
        show = lib[["test", "comparison", "n_called", "n_recovered",
                    "n_wang_in_background", "recall", "fold_enrichment", "pvalue"]]
        print(show.sort_values(["comparison", "test"]).to_string(index=False,
              float_format=lambda x: f"{x:.3g}"))

    ct = df[df["grouping"] != "library"]
    if len(ct):
        print("\n=== cell-type-split, pooled over cell types ===")
        agg = (ct.groupby(["test", "grouping", "comparison"])
                 .agg(n_called=("n_called", "sum"), n_recovered=("n_recovered", "sum"),
                      best_fold=("fold_enrichment", "max"),
                      cell_types=("cell_type", "nunique"))
                 .reset_index())
        print(agg.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

        print("\n=== strongest single cell type per test x comparison ===")
        best = ct.loc[ct.groupby(["test", "grouping", "comparison"])["fold_enrichment"].idxmax()]
        print(best[["test", "grouping", "comparison", "cell_type", "n_called",
                    "n_recovered", "fold_enrichment", "pvalue"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
