# GEI-NJ

GEI-NJ provides helper functions for processing SNV ratio/depth `AnnData` objects,
plotting PCA projections, and constructing neighbor-joining trees from PCA
information.

<div align="center">
    <img src='./share/Schematic_of_GEI-NJ.png' width=1000> 
</div>

## Installation

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

## Contact

* Ruiying Chen (陈睿颖), <chenruiying@genomics.cn>
* Chao Qin (秦超), <qinchao@genomics.cn>