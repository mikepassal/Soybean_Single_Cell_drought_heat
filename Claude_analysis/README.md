# Claude_analysis

Reprocessing, integration and differential expression, answering two questions:

1. **Do we reproduce the collaborator's count matrix?** — yes, bit for bit, wherever the
   FASTQ inputs match. Now confirmed on **all 13 libraries**, not just control.
2. **Is their integration strategy costing us DE signal?** — no. The weak Heat/Drought
   result comes from a **missing replicate-block term** in the DE design.

## Q1: count matrix reproduction — all 13 libraries

`01_compare_count_matrices.py`, our CellRanger **10.1.0** against their **9.0.1**, filtered
matrices, joined on gene *id* and barcode:

| | libraries | result |
|---|---|---|
| rep1 / rep1A / rep3 | 9 | **bit-identical** — 0 differing entries out of 4.4M–13M |
| rep2 | 4 | differ, and the cause is fully diagnosed below |

Nine of thirteen reproduce exactly: same barcodes, same genes, same integer in every cell.
Two CellRanger minor versions apart, that is as strong a confirmation as the comparison can
give. Only one wrinkle, in HD_1: their filtered matrix calls **2 extra barcodes** (~512 UMIs
each, right at the cell-calling threshold) that ours does not. All 7,153 shared cells are
identical, so this is the 9.0.1→10.1.0 cell caller moving a boundary by 0.03%, not a counting
difference.

### The four rep2 libraries are a read-depth artifact on our side, not a discrepancy

Every rep2 library was sequenced **twice** — the 20251117 run and the 20251201 resequencing
run — and the collaborator passed CellRanger both FASTQ paths (hence `_merged` in all four of
their rep2 output directory names). The realignment in `Data/Soybean/realigned_data` passed
only one:

```
--fastqs .../Rep_2/GR0073_20251117_correct_demultiplex_addedRNA/GR0073/untrimmed_MICHAEL_CREATED_for_confirmation
```

Control_2 pins this down arithmetically. The earlier realignment on the Mac used the *other*
half (the 20251201 reseq run) and got 161,279,256 reads; this one used the 20251117 run and
got 164,804,135. They sum to **326,083,391 — exactly** the collaborator's read count for that
library. The two runs are the two halves of one library, and we have been aligning one half at
a time.

The matrices carry the same signature. On shared cells our counts are uniformly *lower*,
never higher:

| library | cells theirs / ours | UMIs retained | 99%+ of differing entries lower in ours | cells with more UMIs in ours |
|---|---|---|---|---|
| Control_2 | 3,310 / 1,742 | 67% | 99.6% | **0** |
| Heat_2 | 26,616 / 20,390 | 74% | 99.6% | **0** |
| Drought_2 | 34,984 / 37,925 | 61% | 99.7% | **0** |
| HD_2 | 5,204 / 4,283 | 78% | 99.3% | **0** |

Per-gene UMI correlation is 0.99993–0.999999 in all four. Not a single cell gains UMIs. A
different aligner, reference or parameter set would scatter differences in both directions;
a missing sequencing run can only subtract. Sequencing saturation confirms it — the rep2
libraries run 41–66% against 74–98% everywhere else, i.e. undersequenced exactly where reads
are missing. The residual 0.3–0.7% of entries that go *up* is ordinary UMI-collapse noise at
lower depth.

Drought_2 is the one place this interacts with cell calling: on ~60% of the reads we call
*more* cells (37,925 vs 34,984), because low saturation flattens the barcode-rank knee. Its
extra 2,970 barcodes are low-UMI (median 310) and 29 of theirs are missing from ours.

**Verdict: we reproduce their count matrix.** The nine libraries with matching inputs match
bit for bit; the four that do not are explained by an input we did not supply, with the
direction, magnitude and read arithmetic all consistent. Re-running rep2 with both `--fastqs`
paths — `realign_fastqs/realign_all_samples.sh` already does this, and was not the script used
for this batch — should close the remaining four.

### What survives into their analysis object

`vs_adata_rna.csv`: of the cells in `adata_rna.h5ad`, we recover **100%** for 11 of 13
libraries, 93.4% for Control_2 and 99.0% for HD_2 — the two rep2 libraries where our shallower
run dropped barcodes below the cell-calling threshold. Their QC is aggressive (e.g. Heat_1:
939 cells kept from 27,618 called), but every cell they kept is a cell we also called.

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

The sample manifest in `pipeline/config.py` detects what is on disk and resolves the data
roots per machine, so the same scripts run unchanged on a partial realignment or on the full
set. All 13 libraries are now present on both sides.

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

**This script has still not been the one that ran.** The batch in
`Data/Soybean/realigned_data` was produced with a single `--fastqs` path pointing at the
20251117 run, so all four rep2 libraries are aligned on roughly half their reads. See Q1
above for the arithmetic. Everything else in that batch is bit-perfect, so rep2 is the only
outstanding item.

## Environment

Two machines have run this. `pipeline/config.py` resolves the data roots by looking for them,
so no edit is needed to move between them.

**Mac (original).** `/Users/michael/Data/Soybean_data`, full collaborator CellRanger tree
(raw + filtered matrices, `metrics_summary.csv`, CellBender h5). Intel, so `annoy` and `torch`
do not build: no BBKNN, Scanorama or scVI locally — Harmony and ComBat cover the two
mechanistically distinct classes. Pins in `.venv`: `numba==0.61.2`, `llvmlite==0.44.0`,
`numpy<2.2`, `pandas<3`, `scanpy==1.11.5`.

**Windows.** Realigned data at `C:/Users/mikep/Data/Soybean/realigned_data`, collaborator
matrices at `C:/Users/mikep/Data/Cellranger_output`, `adata_rna.h5ad` in the repo under
`Data/anndata_export/`. Interpreter is the `Single_cell` conda env
(`C:/Users/mikep/miniconda3/envs/Single_cell/python.exe`; scanpy 1.12.3, anndata 0.13.2,
numpy 2.4.6, pandas 3.0.3). Two caveats:

- The collaborator copy here holds **only the filtered matrices**, flattened to
  `<condition>/<rep>/`. So stage 1 (`metrics_summary.csv`), the raw-matrix comparison and the
  CellBender comparison in script 01 have nothing to run against and are skipped. Script 01
  writes `realigned_metrics_summary.csv` — our side of stage 1 — unconditionally.
- **The env must be activated, not just invoked.** Calling
  `envs/Single_cell/python.exe` directly leaves `envs/Single_cell/Library/bin` off `PATH`,
  MKL's delay-loaded DLLs fail to resolve, and any 2-D `matmul` — so `np.cov`, `np.corrcoef`,
  PCA — kills the interpreter with exit 127 (`0xC06D007F`, delay-load failure) and no
  traceback. `conda activate Single_cell` first, or prepend `Library\bin` to `PATH`. Jupyter
  launched through the activated env is unaffected, which is why the notebooks run fine.
