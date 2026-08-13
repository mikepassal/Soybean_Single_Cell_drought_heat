"""Q2, part 2: run the DE grid over integration variants and grouping schemes.

Consumes the h5ad written by 02_build_and_integrate.py and produces one summary row
per (test x grouping x cell type x contrast), plus the full result tables.

    python Reprocess_and_recluster/scripts/03_run_de.py
    python Reprocess_and_recluster/scripts/03_run_de.py --tests pseudobulk_counts
    python Reprocess_and_recluster/scripts/03_run_de.py --input .../integrated_realigned.h5ad

The headline comparison the summary is built to support:

* ``pseudobulk_counts`` + ``library`` grouping is invariant to every embedding-only
  integration method. If it already shows few Heat/Drought genes, integration is not
  the cause and the answer lies in the counts or the design.
* ``pseudobulk_counts`` split by cell type shows whether the *labels* -- theirs vs
  de-novo, and de-novo from which variant -- change the answer.
* ``pseudobulk_corrected`` and ``wilcoxon_cells`` show what happens when the test is
  fed integration-corrected values instead of raw counts.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as cfg
from pipeline import de as de_mod
from pipeline import integration as integ

warnings.filterwarnings("ignore", category=FutureWarning)

TESTS = ["pseudobulk_counts", "pseudobulk_counts_repblock", "pseudobulk_corrected", "wilcoxon_cells"]


def _corrected_layers(adata: ad.AnnData) -> list[str]:
    """Expression-corrected layers written by the ComBat variants.

    Iterating ``adata.layers`` yields a ``None`` key standing for ``.X``, so the keys
    have to be filtered to real strings before any name test.
    """
    return sorted(l for l in adata.layers if isinstance(l, str) and l.endswith("_corrected"))


def grouping_columns(adata: ad.AnnData, requested: list[str] | None) -> list[str]:
    """Cell-type columns to split by. ``library`` means 'do not split'."""
    candidates = ["library", "celltype_call", "cluster_annot"]
    candidates += [c for c in adata.obs.columns if c.startswith("leiden_")]
    if requested:
        return [c for c in candidates if c in requested]
    return [c for c in candidates if c == "library" or c in adata.obs]


def run_library_level(adata, outdir, test, summaries):
    """One pseudosample per library -- no cell-type split. The reference analysis."""
    print(f"\n### {test} | grouping=library ###", flush=True)

    if test.startswith("pseudobulk_counts"):
        counts, meta = de_mod.pseudobulk(
            adata, ["libraries"], layer="counts", how="sum",
            meta_keys=("treatment", "replicate", "rep_block"),
        )
        counts = de_mod.filter_genes(counts, min_frac_samples=0.5)
        design = "~rep_block + treatment" if test.endswith("repblock") else "~treatment"
        print(f"  {counts.shape[0]} pseudosamples x {counts.shape[1]} genes, design {design}")
        try:
            res = de_mod.deseq(counts, meta, design=design)
        except Exception as exc:                      # rank-deficient design, too few libs
            print(f"  skipped: {exc}")
            return
        tag = f"{test}__library"

    elif test == "pseudobulk_corrected":
        for layer in _corrected_layers(adata):
            vals, meta = de_mod.pseudobulk(
                adata, ["libraries"], layer=layer, how="mean",
                meta_keys=("treatment", "replicate", "rep_block"),
            )
            print(f"  {layer}: {vals.shape[0]} pseudosamples x {vals.shape[1]} genes")
            res = de_mod.linear_model(vals, meta)
            variant = layer.replace("_corrected", "")
            de_mod.save(res, outdir, f"{test}__{variant}__library")
            summaries.append(de_mod.summarize(res, test=test, grouping="library",
                                              variant=variant, cell_type="all"))
            _print(res)
        return

    elif test == "wilcoxon_cells":
        res = de_mod.wilcoxon_cells(adata, layer="lognorm")
        tag = f"{test}__library"

    de_mod.save(res, outdir, tag)
    summaries.append(de_mod.summarize(res, test=test, grouping="library",
                                      variant="-", cell_type="all"))
    _print(res)


def run_celltype_level(adata, outdir, test, group_col, summaries, min_cells, min_reps):
    """Library-level pseudobulk repeated inside each cell type / cluster.

    Pseudosamples stay one-per-library, so replicates remain biological; only the cells
    contributing are restricted. Cell types without enough libraries per treatment are
    skipped rather than modelled on air.
    """
    print(f"\n### {test} | grouping={group_col} ###", flush=True)
    if group_col not in adata.obs:
        return
    sub = adata[adata.obs[group_col].notna()].copy()
    if sub.n_obs == 0:
        return

    layer, how = ("counts", "sum")
    if test == "pseudobulk_corrected":
        corrected = _corrected_layers(adata)
        if not corrected:
            return
        layer, how = corrected[0], "mean"

    counts, meta = de_mod.pseudobulk(
        sub, [group_col, "libraries"], layer=layer, how=how,
        meta_keys=("treatment", "replicate", "rep_block"),
    )
    keep = meta["n_cells"] >= min_cells
    counts, meta = counts.loc[keep.values], meta.loc[keep.values]

    for ct in meta[group_col].unique():
        m = meta[meta[group_col] == ct].copy()
        n_by_t = m["treatment"].value_counts()
        usable = [t for t in cfg.TREATMENTS if n_by_t.get(t, 0) >= min_reps]
        if "Control" not in usable or len(usable) < 2:
            continue
        m = m[m["treatment"].astype(str).isin(usable)]
        c = counts.loc[m.index]
        contrasts = {f"{t} vs Control": ["treatment", t, "Control"] for t in usable if t != "Control"}

        try:
            if test.startswith("pseudobulk_counts"):
                c = de_mod.filter_genes(c, min_frac_samples=0.5)
                design = "~rep_block + treatment" if test.endswith("repblock") else "~treatment"
                if test.endswith("repblock") and m["rep_block"].nunique() < 2:
                    continue
                res = de_mod.deseq(c, m, design=design, contrasts=contrasts)
            else:
                res = de_mod.linear_model(c, m, contrasts=contrasts)
        except Exception as exc:
            print(f"  {ct}: skipped ({type(exc).__name__})")
            continue

        slug = str(ct).replace(" ", "_").replace("/", "-")
        de_mod.save(res, outdir, f"{test}__{group_col}__{slug}")
        summaries.append(de_mod.summarize(res, test=test, grouping=group_col,
                                          variant="-", cell_type=str(ct)))
        n = {k: len(de_mod.sig_genes(v)) for k, v in res.items()}
        print(f"  {str(ct)[:38]:<38} libs={dict(m['treatment'].value_counts()[usable])} sig={n}")


def _print(results):
    for name, res in results.items():
        print(f"    {name:<22} tested={int(res['padj'].notna().sum()):>6,}"
              f"  sig={len(de_mod.sig_genes(res)):>6,}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, default=cfg.RESULTS_ROOT / "integration" / "integrated_collab_h5ad.h5ad")
    ap.add_argument("--tests", nargs="*", default=TESTS, choices=TESTS)
    ap.add_argument("--groupings", nargs="*", default=None)
    ap.add_argument("--min-cells", type=int, default=de_mod.MIN_CELLS_PER_PSEUDOSAMPLE)
    ap.add_argument("--min-reps", type=int, default=de_mod.MIN_REPS_PER_TREATMENT)
    ap.add_argument("--outdir", type=Path, default=cfg.RESULTS_ROOT / "de")
    args = ap.parse_args()

    if not args.input.exists():
        raise SystemExit(f"missing {args.input}\nRun 02_build_and_integrate.py first.")

    adata = sc.read_h5ad(args.input)
    print(adata)
    tables = args.outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    groupings = grouping_columns(adata, args.groupings)
    print(f"\ngroupings: {groupings}")
    print(f"tests:     {args.tests}")

    summaries: list[pd.DataFrame] = []
    for test in args.tests:
        for g in groupings:
            if g == "library":
                run_library_level(adata, tables, test, summaries)
            elif test != "wilcoxon_cells":   # cell-level wilcoxon per cluster is not the comparison here
                run_celltype_level(adata, tables, test, g, summaries, args.min_cells, args.min_reps)

    if not summaries:
        print("\nno results")
        return

    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(args.outdir / "de_summary.csv", index=False)
    print(f"\nwrote {args.outdir / 'de_summary.csv'}  ({len(summary)} rows)")

    print("\n=== library-level (the reference analysis) ===")
    lib = summary[summary["grouping"] == "library"]
    if len(lib):
        print(lib.pivot_table(index=["test", "variant"], columns="comparison",
                              values="n_sig", aggfunc="sum").fillna(0).astype(int).to_string())

    print("\n=== genes called, summed over cell types ===")
    ct = summary[summary["grouping"] != "library"]
    if len(ct):
        piv = ct.pivot_table(index=["test", "grouping"], columns="comparison",
                             values="n_sig", aggfunc="sum").fillna(0).astype(int)
        print(piv.to_string())
        print("\n=== cell types with any called gene ===")
        piv2 = (ct[ct["n_sig"] > 0].pivot_table(index=["test", "grouping"], columns="comparison",
                                                values="cell_type", aggfunc="nunique")
                .fillna(0).astype(int))
        print(piv2.to_string())


if __name__ == "__main__":
    main()
