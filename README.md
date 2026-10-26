# GEI-NJ

GEI-NJ provides an integrated Python workflow for cell relationship analysis from matched SNV ratio and sequencing-depth `AnnData` matrices.

Version **0.2.0.5** uses SNV filtering, binary ALT representation, zero-filling, Manhattan L1 distances, RapidNJ tree construction, and weighted least squares (WLS) rooting. Continuous ALT ratios are also supported. `view_tree_coords()` links the rooted tree to user-defined spatial, PCA, or UMAP coordinates in an offline interactive HTML file.

**Note:** GEI-NJ does not perform SNV calling. Use [GATK Haplotype Caller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), [LineMut Call](https://github.com/RuiyingChenBioinfo/LineMut_v0.2.0) or other suitable workflows to generate the SNV ratio and depth inputs.

## Requirements and installation

Python 3.9+, with `anndata`, `biopython`, `matplotlib`, `numpy`, `pandas`, and `scipy`. These Python dependencies are installed automatically. The coordinate viewer uses bundled HTML and JavaScript; it requires no extra Python package or online JavaScript service. To create a conenvironment and install GEI-NJ:

```bash
conda create -n test_geinj -c conda-forge python=3.12 pip git make cxx-compiler -y
conda activate test_geinj

git clone https://github.com/RuiyingChenBioinfo/GEI-NJ.git
cd GEI-NJ

python -m pip install .

```

RapidNJ is a separate executable and must be on `PATH`, or supplied through `rapidnj_binary`. The benchmark and continuous integration (CI) use [RapidNJ](https://github.com/somme89/rapidNJ) at the following pinned commit. On Linux or macOS with `git`, `make`, and a C++ compiler:

```bash
git clone https://github.com/somme89/rapidNJ.git
git -C rapidNJ checkout ed2d36e219d9db16778b941b5054c0fd021b528a
make -C rapidNJ -j2
export PATH="$PWD/rapidNJ/bin:$PATH"
```

## Input for GEI-NJ

Rows represent cells or cell groups; columns represent SNVs. Depth and ratio objects must have matching, unique observation and SNV names. `comb_depth_ratio()` aligns their order using the ratio object.

- Depths must be finite, nonnegative integer counts.
- At positive depth, ratios must be finite values in `[0, 1]`, calculated as ALT support divided by total support.
- At zero depth, encode the ratio as `0` or `NaN`. Both represent an unobserved entry. Input NaNs produce a warning; NaNs at positive depth raise an error.
- Use original observed ratios, not genotype calls or coarsely rounded ratios. Germline evidence uses exact all-ALT observations (`ratio == 1`).

> [!CAUTION]
> **Invalid depth values, incompatible names, and invalid ratios raise errors. A missing entry does not silently delete an entire SNV. During state assignment, zero-depth entries are marked as NaN in the temporary `state` layer until zero-filling. `X` retains the original ratios.**

## Basic workflow

```python
import anndata as ad
from gei_nj import (
    comb_depth_ratio,
    filter_snvs,
    assign_state,
    zero_filling,
    build_nj,
    view_tree_coords,
)

ad_depth = ad.read_h5ad("unit_depth.h5ad")
ad_ratio = ad.read_h5ad("unit_ratio.h5ad")
adata = comb_depth_ratio(ad_depth, ad_ratio)

adata = filter_snvs(
    adata,
    min_total_depth=5,
    remove_germline=True,
    min_callable_cells_germline=3,
    min_pooled_ratio_germline=0.90,
    min_E_germline=10,
)
adata = assign_state(adata, method="binary", state_low=0.01)
adata = zero_filling(adata, min_callable_depth=1)
tree = build_nj(
    adata,
    wls_weight_power=2,
    do_plot=True,
    use_branch_lengths=False,
    layout="rectangular",
    savepath="geinj_tree.pdf",
    newick_path="geinj_tree.nwk",
)
```

A complete [example](https://github.com/RuiyingChenBioinfo/Analysis_Scripts_for_GEI-NJ_Manuscript/blob/main/02.Nodule_analysis/GEI-NJ_analysis_on_nodule_data.ipynb) and the [relevant input data](https://zenodo.org/records/23665786) are provided. Function signatures and parameter descriptions are provided in the [API reference](./assets/API.html).

### SNV filtering

`filter_snvs()` first excludes SNVs whose depth summed across cells is below `min_total_depth`. With `remove_germline=True`, it also excludes candidate homozygous ALT germline variants satisfying:

```python
remove = (n_callable >= min_callable_cells_germline) & (
    (pooled_ratio >= min_pooled_ratio_germline) | (E >= min_E_germline)
)
```

The three parameters ending in `_germline` apply only when `remove_germline=True`. Here, `n_callable` counts cells with positive depth. `pooled_ratio` is the depth-weighted raw ALT ratio. The e-value combines depth and all-ALT observations across cells; it is an evidence statistic, not a germline probability. This rule identifies candidates with strong all-ALT signals and does not establish their biological origin. It is not intended to detect all heterozygous germline variants. The later `min_callable_depth` setting does not change this filtering rule.

To exclude user-defined SNVs, pass `exclude_snvs` a single-column TSV with **no header**, one complete SNV identifier per line. Identifiers must exactly match `adata.var_names`; chromosome, position, REF, and ALT are not separate columns. For example, `remove_snvs.tsv` contains:

```text
NC_080041.1:122054955:A:G
NC_080041.1:122054973:G:A
NC_080045.1:25890808:G:A
```

Use the identifiers as encoded in your own input matrices. Matching SNVs are removed; entries absent from `adata.var_names` are reported as unmatched. Supply the file with:

```python
adata = filter_snvs(
    adata,
    exclude_snvs="remove_snvs.tsv",
)
```

The supplied [TSV example](./examples/organelle_SNV_candidates.tsv) contains six user-specified loci. It is applied only when explicitly provided and should be replaced with a list appropriate for the dataset. GEI-NJ does not identify or correct RNA editing; known editing sites or other unwanted loci can be excluded this way. Set `remove_germline=False` to disable candidate germline filtering while retaining depth and user-list filtering.

### States and zero-filling

`method="binary"` assigns observed ratios below `state_low=0.01` to 0 and ratios at or above it to 1. `method="continuous"` assigns the raw ALT ratios as continuous states. Both methods use zero-filling at entries with depth below `min_callable_depth`; at the default of 1, only zero-depth entries are filled. Other entries retain their selected state. Zero-filling does not use neighboring cells.

For the continuous alternative, use the same filtered input and rerun state assignment and zero-filling:

```python
adata_cont = assign_state(adata, method="continuous")
adata_cont = zero_filling(adata_cont, min_callable_depth=1)
tree_cont = build_nj(adata_cont, do_plot=False, newick_path="continuous_tree.nwk")
```

`assign_state()` reads the original ratios in `adata.X` and writes a temporary `state` layer. `zero_filling()` writes `zero_filled_state` and `callable`, then removes the temporary state layer. Neither step overwrites `X`. To change the state representation or `min_callable_depth`, rerun `assign_state()` and then `zero_filling()` from the preserved original ratios. Reassigning states clears the previous filled matrix and callability mask; filtering again also clears processed states. Preprocessing functions return a copy by default; `copy=False` updates the object in place.

### Tree construction and plotting

`build_nj()` reads `adata.layers["zero_filled_state"]`, calculates ordinary summed L1 distances between real cells, and calls RapidNJ. It then places a synthetic all-reference profile on this fixed tree by constrained WLS, using weights `(1 + distance_to_reference) ** -wls_weight_power`. The default power is 2; power 0 gives equal weights. The fitted attachment point defines the root, and the returned tree contains only real cells. Negative NJ branch lengths are set to zero before WLS fitting, without collapsing their splits.

Use `layout` to select the tree layout (case-insensitive):

| `layout` | Display |
| --- | --- |
| `"rectangular"` (default) | Horizontal and vertical branches |
| `"triangular"` | Straight diagonal branches |
| `"circular"` | Radial branches and circular arcs |

```python
tree = build_nj(adata, layout="triangular", savepath="geinj_triangular.pdf")
tree = build_nj(adata, layout="circular", figsize=(10, 10), savepath="geinj_circular.pdf")
```

All three layouts support `use_branch_lengths`, `group_key`, `group_colors`, `show_leaf_labels`, and `leafcolor_onside`. In circular plots, side markers form an outer ring and leaf labels follow their radial positions. The layouts use the existing NumPy and Matplotlib dependencies; no additional plotting package is required.

The default `use_branch_lengths=False` draws a topology view with aligned tips. Set it to `True` to draw according to the fitted branch lengths. Layout and branch-length display settings change only the display; the topology, root, and saved Newick branch lengths are unchanged. Branch lengths are derived from the selected SNV representation and distance model; they are not calibrated mutation counts or elapsed time.

Use `group_key` to color cells by an `adata.obs` column. `build_nj()` saves the resolved tip-name-to-color mapping in `tree.tip_colors`, including when `do_plot=False`. Both the static tree plot and `view_tree_coords()` use this same mapping. Without `group_key`, tips are black. `tree.plot_metadata` records the group legend, layout, and branch-length display setting. `do_plot=False` skips plotting, while still constructing the tree and writing `newick_path` if provided. `return_tree_info=True` returns `(tree, info)` instead of the Biopython `Tree` alone.

### View cell coordinates in selected subtrees

Provide a **three-column, tab-separated file without a header** to `coords_path`:

| Column | Content |
| --- | --- |
| 1 | Tip/cell ID, exactly matching a name in `tree.get_terminals()` |
| 2 | Numeric horizontal coordinate, such as spatial x, PC1, or UMAP1 |
| 3 | Numeric vertical coordinate, such as spatial y, PC2, or UMAP2 |

For example, `cell_coordinates.tsv` contains literal tabs between the fields:

```text
Cell_01	12.5	31.0
Cell_02	18.2	27.4
Cell_03	21.8	35.6
```

Include exactly one row for every tip in the tree. Row order is arbitrary: coordinates are matched by ID. Missing, extra, or duplicate IDs, an incorrect column count, and nonnumeric or nonfinite coordinates raise errors. The [coordinate example](./examples/cell_coordinates.tsv) illustrates the format; replace its IDs and coordinates with your own cells.

```python
tree = build_nj(
    adata,
    group_key="cell_type",
    group_colors={"Type A": "#E64B35", "Type B": "#4DBBD5"},
    show_leaf_labels=True,
    savepath="geinj_tree.pdf",
)
html_path = view_tree_coords(
    tree,
    coords_path="cell_coordinates.tsv",
    savepath="geinj_tree_coords.html",
    x_label="X",
    y_label="Y",
    point_size=6.5,
)
```

The example assumes `adata.obs["cell_type"]` exists; adapt the group names to your data. Colors are supplied only to `build_nj()`. The viewer inherits those colors and the selected tree layout and branch-length display setting. For PCA or UMAP, supply the corresponding coordinate file and axis labels:

```python
view_tree_coords(tree, "pca_coordinates.tsv", "tree_pca.html", x_label="PC1", y_label="PC2")
view_tree_coords(tree, "umap_coordinates.tsv", "tree_umap.html", x_label="UMAP1", y_label="UMAP2")
```

`view_tree_coords()` writes the HTML file and returns its `pathlib.Path`. Open the file in a browser. All tips are initially highlighted. Click a tree node to highlight all descendant tips, or a tip/coordinate point to select one cell. Shift-click or enable **Add to selection** to combine selections. Highlighted points retain their tree tip colors; other points are `#e0e0e0`, shown by the **Unselected** legend entry. **Flip X**, **Flip Y**, and **Swap X/Y** adjust the coordinate view; swapping also exchanges the axis labels. The interface and figure labels are in English. X/Y axis lines, tick labels, and axis labels are black; grid lines remain gray.

**Hide tip names** beside **Use branch lengths** in **Rooted tree** hides the tip labels when checked. Tip names are shown by default.

Set `point_size` to a positive finite number to choose the initial coordinate-point radius in SVG plot units (default: `6.5`). The **Point size** control in **Tip coordinates** adjusts all coordinate points, including unselected cells.

**Save current HTML** downloads a standalone file containing the data, colors, current selection, tip-name visibility, point size, coordinate orientation, and interactive controls. The saved file works offline and can be reopened for further selection. Use the tree returned by `build_nj()` for the initial call: Newick stores topology and branch lengths but does not preserve `tip_colors` or `plot_metadata`.

### Processed AnnData structure

| Location | Content |
| --- | --- |
| `adata.X` | Original observed ratios at retained sites, preserving zero-depth 0 or NaN encoding |
| `adata.layers["total_depth"]` | Total sequencing depth for each cell–SNV entry |
| `adata.layers["state"]` | Assigned states; temporary layer removed after zero-filling |
| `adata.layers["callable"]` | Binary mask defined by `min_callable_depth` |
| `adata.layers["zero_filled_state"]` | Zero-filled state matrix used for tree construction |

When custom `layer_name`, `callable_layer`, or `output_layer` names are supplied, the workflow records them in `adata.uns` for the following steps. `build_nj()` uses the recorded filled output layer and does not fall back to `X`.

## Tests

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## Citations

!!!DOI here!!!

## Contact

* Ruiying Chen (陈睿颖), <chenruiying@genomics.cn>
* Chao Qin (秦超), <qinchao@genomics.cn>
* Hai-Xi Sun (孙海汐), <sunhaixi@genomics.cn>
