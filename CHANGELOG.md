# Changelog

## 0.2.0.5

- Set the default binary `state_low` to 0.01 and update the documentation and package version.

## 0.2.0.4.3

- Set the default binary `state_low` to 0.1.
- Add **Hide tip names** beside **Use branch lengths** in the HTML viewer and retain this setting when saving HTML. Tip names remain visible by default.
- Change the viewer accent color to `#45687f`.

## 0.2.0.4.2

- Rename the state method `two-state` to `binary` and set the default `state_low` to 0.01.
- Store the final tip colors on the returned tree, including when static plotting is disabled.
- Add `view_tree_coords()` with a headerless cell-ID/x/y TSV input and inherited tree colors and display settings.
- Link tree-node selection to descendant-tip coordinates in a standalone English HTML viewer. Support multiple selections, axis flips and swapping, and saving the current interactive view.
- Use `#e0e0e0` for unselected cells and add `point_size` plus a saved HTML point-size control. Draw coordinate axes and labels in black while retaining gray grid lines.
- Document user-defined SNV removal lists and coordinate files, and include the HTML template in source and wheel distributions.

## 0.2.0.4.1

- Add the `layout` parameter to `build_nj()`, with rectangular (default), triangular, and circular displays.
- Retain group colors, leaf labels, side markers, and branch-length views in each layout using existing NumPy and Matplotlib dependencies.
- Preserve SNV processing, tree construction, WLS rooting, and Newick output.

## 0.2.0.4

- Use SNV filtering, two-state or continuous ALT ratios, zero-filling, L1 distances, RapidNJ, and WLS rooting.
- Add candidate homozygous ALT germline filtering and optional exclusion by a headerless TSV locus list. Use the `_germline` suffix for its three evidence thresholds.
- Validate matched inputs, depth counts, ratio ranges, and missing-value encoding. Preserve missing zero-depth entries until zero-filling.
- Use `assign_state()` for the default two-state representation or continuous ALT ratio states.
- Keep original ratios in `X`, total depths in `total_depth`, and the tree input in `zero_filled_state`. Remove the temporary `state` layer after zero-filling.
- Place the synthetic all-reference profile on the fixed real-cell NJ tree using constrained WLS.
- Draw aligned-tip topology views by default; retain branch-length views and Newick export.
- Replace the former PCA-based workflow and its processing and tree-building functions.
- Align package versions, add unit tests and CI, and document release archiving.
