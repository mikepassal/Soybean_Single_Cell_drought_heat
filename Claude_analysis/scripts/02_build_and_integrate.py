"""Q2, part 1: build the analysis object and run every integration variant.

Writes one h5ad carrying, side by side, each variant's embedding and its de-novo
leiden clustering, plus a diagnostics table saying how much batch each variant removed
and how much treatment separation survived.

    # all four conditions, from the collaborator's counts (the only complete source now)
    python Reprocess_and_recluster/scripts/02_build_and_integrate.py

    # once heat/drought/HD finish realigning, the same analysis on our own counts
    python Reprocess_and_recluster/scripts/02_build_and_integrate.py --source realigned
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import build, config as cfg
from pipeline import integration as integ


def _sanitize_for_h5ad(adata) -> None:
    """Let anndata write the pandas nullable-string arrays the Seurat export carries.

    They turn up in obs, var and inside ``uns['pca_loadings']``, so chasing them frame
    by frame is fragile. ``False`` (rather than ``True``) writes the pre-0.11
    non-nullable format, which stays readable by older anndata.
    """
    ad.settings.allow_write_nullable_strings = False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default="collab_h5ad", choices=["collab_h5ad", "realigned"])
    ap.add_argument("--variants", nargs="*", default=None, help="subset of variant names")
    ap.add_argument("--resolution", type=float, default=0.8)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--max-pct-chloroplast", type=float, default=None,
                    help="drop cells above this chloroplast %% (default: keep all)")
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "integration")
    ap.add_argument("--tag", default=None, help="suffix for the output filename")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = args.tag or args.source

    print(f"=== loading ({args.source}) ===", flush=True)
    if args.source == "realigned":
        print(cfg.describe_availability("realigned"))
    adata = build.load(args.source)
    print(adata)

    print("\n=== QC ===", flush=True)
    adata = build.qc(adata, max_pct_chloroplast=args.max_pct_chloroplast)
    adata = build.normalize(adata, n_top_genes=args.n_hvg)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes,"
          f" {int(adata.var['highly_variable'].sum())} HVGs")
    print("\ncells per treatment x library:")
    print(pd.crosstab(adata.obs["libraries"], adata.obs["treatment"]).to_string())

    variants = [v for v in integ.VARIANTS if not args.variants or v.name in args.variants]
    # 'theirs' only exists when the collaborator's PCA came along for the ride.
    if args.source != "collab_h5ad":
        variants = [v for v in variants if v.method != "theirs"]

    print("\n=== integration ===", flush=True)
    for v in variants:
        integ.run_variant(adata, v, n_pcs=args.n_pcs, resolution=args.resolution)

    print("\n=== diagnostics ===", flush=True)
    diag = integ.variant_diagnostics(adata, variants)
    print(diag.round(3).to_string(index=False))
    diag.to_csv(args.outdir / f"integration_diagnostics_{tag}.csv", index=False)

    out = args.outdir / f"integrated_{tag}.h5ad"
    # obsp graphs are large and recomputable; drop them before writing.
    for key in list(adata.obsp):
        del adata.obsp[key]
    _sanitize_for_h5ad(adata)
    adata.write_h5ad(out, compression="gzip")
    print(f"\nwrote {out}")
    print(f"wrote {args.outdir / f'integration_diagnostics_{tag}.csv'}")

    print("\nInterpretation guide:")
    print("  batch_purity_repblock -> lower is better mixing of replicate batches")
    print("  treatment_purity      -> HIGHER means condition separation survived")
    print("  a variant that lowers BOTH toward their null values is over-correcting")


if __name__ == "__main__":
    main()
