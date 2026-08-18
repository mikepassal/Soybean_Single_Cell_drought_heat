"""Paths, the sample manifest, and the conventions shared by every script here.

The manifest is the single place that ties together the three views of each library:

* ``run_id``      -- the CellRanger ``--id`` used by the realignment scripts in
                     ``Cluster_analysis_files/Realign_fastqs_cellranger``
* ``collab_dir``  -- the collaborator's CellRanger output directory
* ``library``     -- the ``libraries`` / ``sample`` value inside ``adata_rna.h5ad``

``available_samples()`` looks on disk and returns whatever is actually present, so
the same scripts run on a partial realignment and unchanged once every condition
lands.

Paths are resolved per machine rather than hard-coded, because the analysis has
run on two: the original Mac, and a Windows box whose copy of the collaborator
data holds only the *filtered* matrices, flattened to one directory per library
(``Cellranger_output/<condition>/<rep>/``) with no ``outs/`` wrapper, no raw
matrix, no ``metrics_summary.csv`` and no CellBender h5. ``Sample.matrix_dir``
hides that difference; the stages that need the missing files check for them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# --------------------------------------------------------------------------- paths


def _first_existing(*candidates: Path, default: Path | None = None) -> Path:
    for c in candidates:
        if c.exists():
            return c
    return default if default is not None else candidates[0]


REALIGNED_ROOT = _first_existing(
    Path("/Users/michael/Data/Soybean_data/Realigned_data"),
    Path("C:/Users/mikep/Data/Soybean/realigned_data"),
)
COLLAB_ROOT = _first_existing(
    Path("/Users/michael/Data/Soybean_data/Cell_ranger_output"),
    Path("C:/Users/mikep/Data/Cellranger_output"),
)
COLLAB_H5AD = _first_existing(
    Path("/Users/michael/Data/Soybean_data/adata_rna.h5ad"),
    REPO_ROOT / "Data" / "anndata_export" / "adata_rna.h5ad",
)

# The Windows copy is flattened: <condition>/<rep>/matrix.mtx.gz, filtered only.
# The Mac copy is the full CellRanger tree: <long_condition_dir>/<run_dir>/outs/.
COLLAB_LAYOUT = "flat" if (COLLAB_ROOT / "control" / "rep1" / "matrix.mtx.gz").exists() else "cellranger"

ANALYSIS_ROOT = REPO_ROOT / "Claude_analysis"
RESULTS_ROOT = ANALYSIS_ROOT / "results"                 # gitignored: ~3 GB of tables
SUMMARY_ROOT = ANALYSIS_ROOT / "results_summary"         # committed: the small tables
CACHE_ROOT = RESULTS_ROOT / "cache"

# Collaborator CellRanger runs are grouped into one parent directory per condition.
COLLAB_CONDITION_DIRS = {
    "control": "cellranger_out_n_Cellbender_out_controls_C1_C1A_C2_control3_andrewgenome",
    "heat": "cellranger_out_n_Cellbender_out_heat_H1_H2_heat3_andrewgenome_2",
    "drought": "cellranger_out_n_Cellbender_out_drought_D1_D2_drought3_andrewgenome",
    "heat_drought": "cellranger_out_n_Cellbender_out_HD_HD1_HD2_HD3_andrewgenome",
}

# CellBender was run per condition into its own subdirectory of the above.
COLLAB_CELLBENDER_DIRS = {
    "control": "cellbender_auto_Control_v2",
    "heat": "cellbender_auto_heat_andygenome_v3",
    "drought": "cellbender_auto_drought_v3",
    "heat_drought": "cellbender_auto_HD_v3_andrewgenome",
}

# Condition directory names under the flat layout.
COLLAB_FLAT_CONDITION_DIRS = {
    "control": "control",
    "heat": "heat",
    "drought": "drought",
    "heat_drought": "HD",
}

# --------------------------------------------------------------------------- labels

# Treatment labels as they appear in adata_rna.h5ad. Control first so it stays the
# DESeq2 reference level everywhere.
TREATMENTS = ["Control", "Heat", "Drought", "HD"]
CONDITION_TO_TREATMENT = {
    "control": "Control",
    "heat": "Heat",
    "drought": "Drought",
    "heat_drought": "HD",
}
CONTRASTS = {f"{t} vs Control": ["treatment", t, "Control"] for t in TREATMENTS[1:]}

# Organelle contigs in the soy_Williams82_v3_with_organelles reference. These are
# absent from adata_rna.h5ad (the collaborator dropped them after computing pctCP)
# but present in every CellRanger matrix, and chloroplast load differs sharply
# between libraries, so they get their own QC treatment.
CHLOROPLAST_PREFIX = "GlmaCp"
MITO_PREFIXES = ("GlmaCt", "GlmaCr")


@dataclass(frozen=True)
class Sample:
    run_id: str  # CellRanger --id in the realignment
    condition: str  # control / heat / drought / heat_drought
    replicate: str  # rep1 / rep1A / rep2 / rep3
    collab_dir: str  # CellRanger output dir name in the collaborator tree
    library: str  # `libraries` value in adata_rna.h5ad

    @property
    def treatment(self) -> str:
        return CONDITION_TO_TREATMENT[self.condition]

    @property
    def realigned_outs(self) -> Path:
        return REALIGNED_ROOT / self.condition / self.run_id / "outs"

    @property
    def collab_outs(self) -> Path:
        return COLLAB_ROOT / COLLAB_CONDITION_DIRS[self.condition] / self.collab_dir / "outs"

    @property
    def collab_flat_dir(self) -> Path:
        """The one directory holding this library's collaborator matrix, flat layout."""
        return COLLAB_ROOT / COLLAB_FLAT_CONDITION_DIRS[self.condition] / self.replicate.lower()

    @property
    def cellbender_h5(self) -> Path:
        base = COLLAB_ROOT / COLLAB_CONDITION_DIRS[self.condition]
        cb = base / COLLAB_CELLBENDER_DIRS[self.condition] / self.collab_dir
        return cb / f"cellbender_{self.collab_dir}_filtered.h5"

    def matrix_dir(self, source: str, filtered: bool = True) -> Path:
        """``source`` is 'realigned' or 'collaborator'.

        Under the flat collaborator layout only the filtered matrix exists; the
        returned raw path will simply not be there, which callers already test for.
        """
        if source == "collaborator" and COLLAB_LAYOUT == "flat":
            return self.collab_flat_dir if filtered else self.collab_flat_dir / "raw_feature_bc_matrix"
        outs = self.realigned_outs if source == "realigned" else self.collab_outs
        return outs / ("filtered_feature_bc_matrix" if filtered else "raw_feature_bc_matrix")

    def has(self, source: str) -> bool:
        return (self.matrix_dir(source) / "matrix.mtx.gz").exists()


SAMPLES: list[Sample] = [
    # control ---------------------------------------------------------------
    Sample("Control_1", "control", "rep1", "RNA_leaf_C1_GEX", "leaf_control_rep1"),
    Sample("Control_1A", "control", "rep1A", "RNA_leaf_C1A_GEX", "leaf_control_rep1A"),
    Sample("Control_2", "control", "rep2", "RNA_leaf_C2_merged_GEX", "leaf_control_rep2"),
    Sample("Control_3", "control", "rep3", "RNA_leaf_control_3_GEX", "leaf_control_rep3"),
    # heat ------------------------------------------------------------------
    Sample("Heat_1", "heat", "rep1", "RNA_leaf_H1_GEX_andrewgenome", "leaf_Heat_rep1"),
    Sample("Heat_2", "heat", "rep2", "RNA_leaf_H2_merged_GEX_andrewgenome", "leaf_Heat_rep2"),
    Sample("Heat_3", "heat", "rep3", "RNA_leaf_heat_3_GEX_andrewgenome", "leaf_Heat_rep3"),
    # drought ---------------------------------------------------------------
    Sample("Drought_1", "drought", "rep1", "RNA_leaf_D1_GEX_andrewgenome", "leaf_drought_rep1"),
    Sample("Drought_2", "drought", "rep2", "RNA_Leaf_D2_merged_GEX_andrewgenome", "leaf_drought_rep2"),
    Sample("Drought_3", "drought", "rep3", "RNA_leaf_drought_3_GEX_andrewgenome", "leaf_drought_rep3"),
    # heat + drought --------------------------------------------------------
    Sample("HD_1", "heat_drought", "rep1", "RNA_leaf_HD1_GEX_andrewgenome", "leaf_HD_rep1"),
    Sample("HD_2", "heat_drought", "rep2", "RNA_Leaf_HD2_merged_GEX_andrewgenome", "leaf_HD_rep2"),
    Sample("HD_3", "heat_drought", "rep3", "RNA_leaf_HD_3_GEX_andrewgenome", "leaf_HD_rep3"),
]

SAMPLES_BY_RUN_ID = {s.run_id: s for s in SAMPLES}
SAMPLES_BY_LIBRARY = {s.library: s for s in SAMPLES}


def available_samples(source: str = "realigned", samples: list[Sample] | None = None) -> list[Sample]:
    """Samples whose matrix for ``source`` is on disk right now.

    This is what makes the scripts rerunnable: today it returns the four control
    libraries, and it will return all thirteen once the remaining CellRanger runs
    finish, with no edits anywhere.
    """
    return [s for s in (samples or SAMPLES) if s.has(source)]


def missing_samples(source: str = "realigned") -> list[Sample]:
    return [s for s in SAMPLES if not s.has(source)]


def barcode_to_collab_obs_name(barcode: str, library: str) -> str:
    """``AAACAGCCATGTCGCG-1`` -> ``leaf_control_rep1_AAACAGCCATGTCGCG-1``."""
    return f"{library}_{barcode}"


def describe_availability(source: str = "realigned") -> str:
    have, missing = available_samples(source), missing_samples(source)
    lines = [f"{source}: {len(have)}/{len(SAMPLES)} libraries on disk"]
    for cond in COLLAB_CONDITION_DIRS:
        got = [s.run_id for s in have if s.condition == cond]
        lines.append(f"  {cond:<13} {len(got)}/{sum(s.condition == cond for s in SAMPLES)}"
                     + (f"  [{', '.join(got)}]" if got else "  (none yet)"))
    return "\n".join(lines)
