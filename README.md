# GEI-NJ

GEI-NJ provides helper functions for processing SNV ratio/depth `AnnData` objects,
plotting PCA projections, and constructing neighbor-joining trees from PCA
coordinates.

## Installation

Install locally while developing:

```bash
pip install -e .
```

Install directly from GitHub after publishing:

```bash
pip install git+https://github.com/RuiyingChenBioinfo/GEI-NJ.git
```

## Usage

```python
from gei_nj import (
    add_callable_layer,
    build_and_plot_nj,
    comb_depth_ratio,
    filt_adata_by_mutation_nnz,
    filt_snv_by_depth_ratio,
    impute_non_callable_by_snv_mean,
    plot_snv_pca_by_group,
    revise_snv_ratio,
)

adata = comb_depth_ratio(ad_depth, ad_ratio)
adata = filt_snv_by_depth_ratio(adata, min_tot_depth=20)
adata = revise_snv_ratio(adata)
adata = add_callable_layer(adata, min_callable_depth=3)
adata = impute_non_callable_by_snv_mean(adata)

plot_snv_pca_by_group(adata, group_by="group")
tree = build_and_plot_nj(adata, group_key="group")
```

## Modules

- `gei_nj.process_data`: SNV depth/ratio filtering, discretization, callability,
  and imputation helpers.
- `gei_nj.construct_nj`: PCA plotting and neighbor-joining tree construction.

## Publishing To GitHub

From this folder:

```bash
git init
git add .
git commit -m "Package GEI-NJ utilities"
git branch -M main
git remote add origin https://github.com/RuiyingChenBioinfo/GEI-NJ.git
git push -u origin main
```