import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import scanpy as sc 
    import matplotlib.pyplot as plt
    import seaborn as sns

    return (sc,)


@app.cell
def _(sc):
    orig_rna_soybean= sc.read_h5ad(r"C:\Users\mikep\git\Soybean_Single_Cell_drought_heat\Data\anndata_export\adata_rna_scaled.h5ad")
    integrated_soybean = sc.read_h5ad(r"C:\Users\mikep\git\Soybean_Single_Cell_drought_heat\Data\anndata_export\adata_integrated.h5ad")

    return integrated_soybean, orig_rna_soybean


@app.cell
def _(integrated_soybean, orig_rna_soybean):
    all(orig_rna_soybean.obs.index == integrated_soybean.obs.index)
    return


@app.cell
def _(orig_rna_soybean):
    orig_rna_soybean.obs
    return


@app.cell
def _(orig_rna_soybean):
    orig_rna_soybean.var.head()
    return


@app.cell
def _(integrated_soybean):
    integrated_soybean.obs
    return


if __name__ == "__main__":
    app.run()
