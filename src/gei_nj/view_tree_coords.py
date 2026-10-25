def view_tree_coords(
    tree,
    coords_path,
    savepath="geinj_tree_coords.html",
    x_label="X",
    y_label="Y",
    point_size=6.5,
):
    """Link a rooted tree to cell coordinates in a standalone HTML file.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Tree returned by build_nj, including its saved tip_colors.
    coords_path : str or Path
        Three-column, headerless TSV: tip ID, x, y. IDs must match the
        tree tips exactly; row order does not matter.
    savepath : str or Path
        Output HTML. No internet connection is needed to open the file.
    x_label, y_label : str
        Axis labels, for example X/Y, PC1/PC2 or UMAP1/UMAP2.
    point_size : float
        Initial coordinate-point radius in plot units. Adjustable in HTML.

    Returns
    -------
    pathlib.Path
        Saved HTML path. The viewer inherits tree colors and layout.

    Notes
    -----
    Click a node to highlight its descendant tips. Shift-click or enable
    Add to selection for a union of selections. Flip X, Flip Y and Swap X/Y
    affect only the displayed coordinates. Save current HTML preserves
    the data, selection, colors, point size, tip-name visibility and axis settings.
    """
    import json
    from pathlib import Path
    import numpy as np
    import pandas as pd

    if not np.isfinite(point_size) or point_size <= 0:
        raise ValueError("point_size must be positive and finite.")
    if not hasattr(tree, "tip_colors"):
        raise ValueError("tree has no saved tip colors; use the tree returned by build_nj.")
    table = pd.read_csv(coords_path, sep="\t", header=None, dtype=str,
                        keep_default_na=False)
    if table.shape[1] != 3:
        raise ValueError("coords_path must contain three headerless columns: tip ID, x, y.")
    names = [tip.name for tip in tree.get_terminals()]
    if table.iloc[:, 0].duplicated().any():
        raise ValueError("Coordinate tip IDs must be unique.")
    ids = set(table.iloc[:, 0])
    missing, extra = set(names) - ids, ids - set(names)
    if missing or extra:
        raise ValueError(f"Coordinate IDs must match tree tips: {len(missing)} missing, {len(extra)} extra.")
    coordinates = table.iloc[:, 1:].to_numpy(dtype=float)
    if not np.isfinite(coordinates).all():
        raise ValueError("Coordinates must be finite numbers.")
    coords = dict(zip(table.iloc[:, 0], coordinates.tolist()))
    node_ids = {node: f"node{i}" for i, node in enumerate(tree.find_clades())}

    def encode(node):
        item = dict(id=node_ids[node], length=float(node.branch_length or 0.))
        if node.clades:
            item["label"] = str(node.name or ("Root" if node is tree.root else node_ids[node]))
            item["children"] = [encode(child) for child in node.clades]
        else:
            item["tip_id"] = node.name
        return item

    metadata = tree.plot_metadata
    data = dict(tree=encode(tree.root), tip_colors=tree.tip_colors, coords=coords,
                groups=[dict(label=str(group), color=color)
                        for group, color in metadata["group_colors"].items()],
                labels=[str(x_label), str(y_label)], layout=metadata["layout"])
    state = dict(selected=[node_ids[tree.root]], additive=False,
                 branchLengths=bool(metadata["use_branch_lengths"]), hideTipNames=False,
                 flipX=False, flipY=False, swapXY=False,
                 pointSize=float(point_size), views=None)
    template_path = Path(__file__).with_name("templates") / "tree_coords.html"
    html = template_path.read_text(encoding="utf-8")
    html = html.replace("__GEINJ_STATE__", json.dumps(state))
    html = html.replace("__GEINJ_DATA__", json.dumps(data, ensure_ascii=False,
                        allow_nan=False).replace("<", "\\u003c"))
    path = Path(savepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path
