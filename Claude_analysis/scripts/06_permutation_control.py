"""Negative control: does adding `rep_block` manufacture significance?

The worry is reasonable. The 13 libraries are simultaneously the pseudosamples (the
replicates DESeq2 fits) and the things `rep_block` is derived from, so it can look like
the same data is being spent twice.

It is not circular -- `rep_block` comes from FASTQ provenance, which is fixed before any
expression is examined, and each library still contributes exactly one observation
carrying one treatment label and one block label. But "not circular in principle" is not
proof, so this script measures it.

**Restricted permutation.** Treatment labels are shuffled *within* each block. That
destroys the treatment effect while preserving the block structure exactly, which is the
correct null for a randomized complete block design -- an unrestricted shuffle would also
scramble the blocks and test the wrong hypothesis.

If `~rep_block + treatment` were inventing signal, it would call genes on permuted labels
at a similar rate to real labels. If instead it calls ~nothing on permuted labels while
calling 106/110 on the real ones, the extra genes are a power gain, not inflation.

    python Reprocess_and_recluster/scripts/06_permutation_control.py --n-perm 10
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as cfg
from pipeline import de as de_mod

warnings.filterwarnings("ignore")

DESIGNS = ["~treatment", "~rep_block + treatment"]


def permute_within_block(meta: pd.DataFrame, rng: np.random.Generator) -> pd.Series:
    """Shuffle treatment labels among the libraries inside each block."""
    out = meta["treatment"].astype(str).copy()
    for block, idx in meta.groupby("rep_block", observed=True).groups.items():
        idx = list(idx)
        out.loc[idx] = rng.permutation(out.loc[idx].to_numpy())
    return out


def count_sig(counts, meta, design) -> dict[str, int]:
    try:
        res = de_mod.deseq(counts, meta, design=design, n_cpus=4)
    except Exception as exc:
        print(f"    fit failed ({type(exc).__name__}), skipping", flush=True)
        return {}
    return {k: len(de_mod.sig_genes(v)) for k, v in res.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path,
                    default=cfg.RESULTS_ROOT / "integration" / "integrated_collab_h5ad.h5ad")
    ap.add_argument("--n-perm", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "de")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    counts, meta = de_mod.pseudobulk(
        adata, ["libraries"], layer="counts", how="sum",
        meta_keys=("treatment", "replicate", "rep_block"),
    )
    counts = de_mod.filter_genes(counts, min_frac_samples=0.5)
    print(f"{counts.shape[0]} libraries x {counts.shape[1]} genes")
    print(f"\nlibraries per block x treatment:")
    print(pd.crosstab(meta["rep_block"], meta["treatment"]).to_string())

    print("\n=== observed (real labels) ===", flush=True)
    observed = {}
    for design in DESIGNS:
        observed[design] = count_sig(counts, meta, design)
        print(f"  {design:<24} {observed[design]}")

    print(f"\n=== {args.n_perm} within-block permutations ===", flush=True)
    rng = np.random.default_rng(args.seed)
    rows = []
    for i in range(args.n_perm):
        pm = meta.copy()
        pm["treatment"] = pd.Categorical(
            permute_within_block(meta, rng), categories=cfg.TREATMENTS
        )
        # A shuffle can leave a treatment absent from a block; skip those draws.
        if pm.groupby("rep_block", observed=True)["treatment"].nunique().min() < 2:
            continue
        for design in DESIGNS:
            sig = count_sig(counts, pm, design)
            for contrast, n in sig.items():
                rows.append({"perm": i, "design": design, "comparison": contrast, "n_sig": n})
        done = {d: {r["comparison"]: r["n_sig"] for r in rows if r["perm"] == i and r["design"] == d}
                for d in DESIGNS}
        print(f"  perm {i + 1:>2}: {done['~rep_block + treatment']}", flush=True)

    if not rows:
        raise SystemExit("no usable permutations")

    null = pd.DataFrame(rows)
    null.to_csv(args.outdir / "permutation_null.csv", index=False)

    print("\n=== null distribution of significant genes ===")
    summary = (null.groupby(["design", "comparison"])["n_sig"]
                   .agg(["mean", "median", "max", "size"])
                   .rename(columns={"size": "n_perm"}))
    summary["observed"] = [
        observed[d].get(c, np.nan) for d, c in summary.index
    ]
    # Permutation p-value: how often the null reaches the observed count.
    summary["p_perm"] = [
        (null[(null.design == d) & (null.comparison == c)]["n_sig"] >= observed[d].get(c, np.inf)).mean()
        for d, c in summary.index
    ]
    print(summary.round(3).to_string())

    print(f"\nwrote {args.outdir / 'permutation_null.csv'}")
    print("\nRead: if the block term were manufacturing significance, the null means for")
    print("'~rep_block + treatment' would approach its observed counts. They should not.")


if __name__ == "__main__":
    main()
