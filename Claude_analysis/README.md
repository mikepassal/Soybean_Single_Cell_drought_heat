# Claude_analysis

Reprocessing, integration and differential expression, answering two questions:

1. **Do we reproduce the collaborator's count matrix?** — yes, bit for bit, wherever the
   FASTQ inputs match.
2. **Is their integration strategy costing us DE signal?** — no. The weak Heat/Drought
   result comes from a **missing replicate-block term** in the DE design.

## Headline result

`rep_block` (rep1/rep2/rep3 — three complete passes through the 4-condition experiment) is
the dominant axis of variation: it explains **72% of PC1**, which is 46% of all library-level
variance, against treatment's 13%. Adding it to the design:

| Contrast | `~treatment` | `~rep_block + treatment` |
|---|---|---|
| Heat vs Control | 8 | **106** |
| Drought vs Control | 12 | **110** |
| HD vs Control | 2,398 | 4,199 |

Validated against Wang et al. (independent bulk RNA-seq): Heat enrichment *rises* 4.0x → 6.4x
(p 0.22 → 2.7e-10) as genes are added, so the extra calls are real signal, not noise.

Running DE on integration-corrected values instead of raw counts recovers **zero** validated
Heat genes under every grouping tried — the likely explanation for the collaborator's result
if their `FindMarkers` read the integrated assay rather than the RNA assay.

## Layout

```
Claude_analysis/
  pipeline/            importable modules (config, matrices, build, integration, de, reference)
  scripts/             01-06, run in order
  realign_fastqs/      corrected CellRanger array job (both rep2 FASTQ paths + preflight)
  results_summary/     small summary tables -- tracked in git
  results/             full output, ~3 GB -- gitignored
```

## Running it

Everything runs from `.venv/bin/python` at the repo root (see Environment below).

```bash
.venv/bin/python Claude_analysis/scripts/01_compare_count_matrices.py
.venv/bin/python Claude_analysis/scripts/02_build_and_integrate.py
.venv/bin/python Claude_analysis/scripts/03_run_de.py
.venv/bin/python Claude_analysis/scripts/04_benchmark_vs_wang.py
.venv/bin/python Claude_analysis/scripts/05_cluster_composition.py
.venv/bin/python Claude_analysis/scripts/06_permutation_control.py
```

The sample manifest in `pipeline/config.py` detects what is on disk, so these run today on
the four control libraries and unchanged once heat, drought and HD finish realigning.

## Experimental structure

The experiment was run three times over, each pass covering all four conditions. That is
what rep1/rep2/rep3 are, and it is visible in the FASTQ provenance: rep2 is one sequencing
submission (`GR0073`) holding all four treatments, rep3 likewise (`GR0192`). `rep1A` is a
second control library from the rep1 batch with no counterpart elsewhere, so it folds into
rep1.

```
             rep1  rep1A  rep2  rep3
Control      3740   3452   823  1683
Heat          939     --  2381  2193
Drought      1434     --  2855  1036
HD           5538     --  3691   702
```

This is a randomized complete block design. `rep_block` is **crossed** with treatment (every
block holds all four conditions), so both effects are estimable. `libraries` is **nested**
inside treatment (each library is one condition), so it can never be used as a batch
covariate or integration key without removing the biology along with the batch.

## What the collaborator's object actually is

`adata_rna.h5ad` is CellRanger 9.0.1 → CellBender → Seurat 5.3.0 SCTransform → **integration
run separately within each of the four conditions** → merge → annotate. The four
`integration_final.rds` files named in `obs['source_rds']` are the evidence; the merged
object then got a shared 3000-gene PCA and joint clustering, which is where `celltype_call`
and `cluster_annot` come from.

Its `.X` is the RNA log-normalised layer and `layers['counts']` is raw counts
(`uns['source']['X_is']`), so pseudobulk off the counts layer is honest. The SCT and
integrated values are not in the file.

## Design of the comparison

### Integration variants (`pipeline/integration.py`)

Two axes, varied independently:

**What the batch key is.** `libraries` is *nested* inside `treatment` — every cell in a
library shares a condition — so correcting on it asks the method to remove variation
perfectly confounded with the biology. `rep_block` (rep1/rep2/rep3, rep1A folded into rep1)
is *crossed* with treatment: every block contains all four conditions.

```
             Control  Heat  Drought  HD
rep1            3740   939     1434  5538
rep2             823  2381     2855  3691
rep3            1683  2193     1036   702
```

**What the method corrects.** Harmony moves an *embedding* only, so raw-count pseudobulk DE
is mathematically unaffected by it — it can only reach DE by changing cluster assignment.
ComBat rewrites *expression values*, so it propagates into any DE run on corrected values.
This is the distinction that decides whether "integration reduced the signal" is even
mechanically possible for a given pipeline.

| variant | method | batch key | expectation |
|---|---|---|---|
| `unintegrated` | — | — | baseline |
| `harmony_repblock` | embedding | crossed | defensible: batch removed, contrast preserved |
| `harmony_library` | embedding | nested | condition separation should collapse |
| `combat_repblock` | expression | crossed | corrected values, contrast preserved |
| `combat_library` | expression | nested | closest analogue to DE on a Seurat integrated assay |
| `theirs` | — | — | their PCA and `cluster_annot`, as shipped |

Diagnostics report `batch_purity_repblock` (lower = better mixing) alongside
`treatment_purity` (higher = condition separation survived). **A variant that drives both
toward their null values is over-correcting** — that is the failure mode being hunted.

### DE grid (`pipeline/de.py`)

*What the test consumes*

- `pseudobulk_counts` — sum raw counts per library, DESeq2. Replicates are libraries, i.e.
  real biological replicates. The reference analysis.
- `pseudobulk_counts_repblock` — same, design `~ rep_block + treatment`, absorbing the
  replicate batch that chloroplast load tracks.
- `pseudobulk_corrected` — mean ComBat-corrected expression per library, then OLS. DESeq2's
  negative binomial is wrong for corrected values (continuous, can be negative), so this arm
  uses a linear model. This is the analogue of DE on a Seurat integrated assay.
- `wilcoxon_cells` — cell-level rank-sum on log-normalised values. What `FindMarkers` does.
  Included as a comparator, not a recommendation: cells within a library are not
  independent, so the p-values are severely anti-conservative.

*How cells are grouped* — `library` (no split), `celltype_call` and `cluster_annot` (theirs),
and `leiden_<variant>` (de novo, one per integration variant).

### Why Wang et al. is the arbiter (`04_benchmark_vs_wang.py`)

Ranking strategies by "how many genes did it call" just rewards the most anti-conservative
one. `Referenced_data/Wang_DE_gene_lists/` holds independent bulk RNA-seq calls for the same
three stresses at 24h; 66–78% of those genes are present in our data after mapping
Wm82.a4.v1 → Glyma 1.1. Hypergeometric enrichment against them, with the background
restricted to genes actually tested, asks which strategy recovers real biology.

## Known confound: chloroplast load

`pctCP` varies ~20x across libraries and tracks replicate batch, not condition (median 0.31
Control / 0.16 Heat / 0.63 Drought / **5.85 HD**; all rep3 libraries <0.5). HD is dominated
by the high-chloroplast rep1+rep2 libraries. Chloroplast transcripts inflate library size, so
`normalize_total` deflates every nuclear gene — and the existing HD result is 1,383 up /
3,624 down, which is that artifact's signature.

So there are two hypotheses on the table, not one: integration suppressing Heat/Drought, and
**HD being partly inflated**. `--max-pct-chloroplast` and the `~ rep_block + treatment`
design exist to separate them.

Organelle genes (`GlmaCp*`, `GlmaCt*`/`GlmaCr*`) were dropped from `adata_rna.h5ad` after
`pctCP` was computed but are present in the CellRanger matrices, so the `realigned` source
can test this directly once all conditions land.

## Additional checks

### Negative control (`06_permutation_control.py`)

Since the 13 libraries are simultaneously the pseudosamples and the source of `rep_block`,
this shuffles treatment labels **within** each block — destroying the treatment effect while
preserving block structure, the correct null for a blocked design — and checks that the
block term is not manufacturing significance.

Result over 10 permutations (genes called under `~rep_block + treatment`):

| Contrast | null mean | null median | null max | observed |
|---|---|---|---|---|
| Heat vs Control | 4.3 | 1 | 30 | **106** |
| Drought vs Control | 9.0 | 0.5 | 70 | **110** |
| HD vs Control | 3.2 | 0 | 14 | **4,199** |

No permutation reached the observed count. The block term *is* modestly more liberal under
the null than `~treatment` (whose null means are 0.2 / 0.3 / 0.0), so it is not free — but the
observed signal sits 12-25x above the null mean and outside the null range entirely. With only
10 permutations the p-value resolution is coarse (0/10); the separation in effect size is the
stronger evidence.

### Cluster composition (`05_cluster_composition.py`)

Reports how batch-confounded each clustering is. Motivated by a concrete case: their
`Mesophyll-7` is 91% rep3 cells, so each condition contributes only its single rep3 library
and the stress contrasts silently drop out of a per-cell-type DE rather than returning a null.

### Realignment fix (`realign_fastqs/realign_all_samples.sh`)

Replaces the four per-condition scripts in `Cluster_analysis_files/`. Those passed only one
`--fastqs` path for the rep2 libraries; the collaborator passed two (the 20251117 run and the
20251201 resequencing run — hence `_merged` in their output directory names). Includes a
preflight check that aborts loudly rather than silently half-aligning a library.

Note: the rep2 mismatch is **expected** in our hands — the first run sat behind a symlink that
was not accessible in the original download.

## Environment

Intel Mac, so `annoy` and `torch` do not build: no BBKNN, Scanorama or scVI locally. Harmony
and ComBat cover the two mechanistically distinct classes. Pins in
`pipeline/../../.venv`: `numba==0.61.2`, `llvmlite==0.44.0`, `numpy<2.2`, `pandas<3`,
`scanpy==1.11.5`.
