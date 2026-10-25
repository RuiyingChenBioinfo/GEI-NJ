# Validation history of GEI-NJ

## Binary default update 0.2.0.5

- Round 1: set the binary `state_low` default to 0.01 and update current version references and documentation. All 49 pytest cases passed with real RapidNJ, including dense/sparse boundaries at 0.009 and 0.01, explicit threshold overrides, missing values, callability, and continuous-state restoration.
- Round 2: independently reviewed the source, tests, and documentation, installed the rebuilt wheel, and passed all 49 cases again. Confirmed version 0.2.0.5, the installed default of 0.01, and matching source/installed Python files and HTML template. Tree construction, the coordinate viewer, and dependencies are unchanged. Earlier version references below and in migration notes describe historical releases.

## Binary default adjustment within 0.2.0.4.3

- Round 1: changed only the `assign_state()` default from `state_low=0.01` to `0.1` in production code and updated current documentation. All 49 pytest cases passed with real RapidNJ, including dense/sparse boundaries at 0.099 and 0.1, explicit threshold overrides, missing values, callability, and continuous-state restoration.
- Round 2: independently reviewed source and documentation changes, installed the rebuilt wheel, and passed all 49 cases again. Confirmed the installed default is 0.1 and all installed Python files and the HTML template match the source. Version remains 0.2.0.4.3; viewer, construction code, and dependencies are unchanged.

## Viewer update 0.2.0.4.3

- Round 1: reviewed the minimal viewer/state changes and ran all 49 pytest cases successfully with real RapidNJ. Version metadata, documentation, and the version-alignment test were updated together.
- Checked HTML from a 32-cell, 123-SNV RapidNJ/WLS fixture in Chromium for rectangular, triangular, and circular layouts. Tip names start visible; hiding them leaves node positions, markers, coordinate points, and selection unchanged. Node selection and branch-length toggling continue to work.
- Saved and reopened each layout offline with tip names hidden. Visibility, selection, point size, and coordinate transforms were retained; names could be shown again and coordinate points remained interactive. No browser script errors occurred.
- Verified `#45687f` for the save button (including hover), selected internal nodes, selection text, and checkbox accents. Reviewed the circular layout and confirmed that the additional control fits the mobile page width.
- Round 2: independently reviewed the source changes, installed the built wheel, and passed all 49 pytest cases again. Installed Python files and the HTML template match the reviewed source. Preprocessing, NJ/WLS construction, dependencies, and function signatures are unchanged.

## Viewer refinements within 0.2.0.4.2

- Round 1: renamed the public function, module, imports, examples, and tests to `view_tree_coords`; added the initial `point_size` parameter and HTML control. All 49 pytest cases passed with the pinned RapidNJ executable. Independently reviewed the code and documentation changes.
- Round 2: installed the built wheel and passed all 49 cases again. Confirmed that the installed package exports the new function, omits the old module, and includes the current HTML template. The preprocessing and NJ/WLS source files, dependencies, and packaging configuration are unchanged.
- Tested HTML generated from a 32-cell, 123-SNV RapidNJ/WLS fixture in Chromium, both from source and from the installed wheel. Checked `#e0e0e0` for unselected cells and the legend, black axes/ticks/labels, and unchanged gray grid lines.
- Changed point size in the viewer and verified that all coordinate points resized while tree nodes, coordinates, colors, and selection were preserved. Saved and reopened the HTML offline; point size, selection, and axis transforms were restored, and further resizing and selection worked without script errors.

## Tree-coordinate update 0.2.0.4.2

Two local review rounds were completed for this update.

- Round 1: reviewed the default threshold/method changes, shared tip-color mapping, coordinate parser, HTML viewer, and documentation. All 48 pytest cases passed with the real pinned RapidNJ executable. The NJ, WLS-rooting, and existing tree-coordinate calculation functions are unchanged. Static leaf-label colors matched `tree.tip_colors` in all three layouts, including a semitransparent custom color.
- Round 2: built the source distribution and wheel, installed the wheel in place of the editable package, and passed the same 48 cases. Confirmed that imports and template loading used the installed package. The HTML template is present in both distributions; all six Python runtime dependencies remain declared.
- In Chromium, generated HTML from a 32-cell RapidNJ/WLS fixture and tested rectangular, triangular, and circular layouts. Checked node descendants, Shift/additive selection, point-to-tip selection, exact inherited colors, and `#F0F0F0` for unselected points.
- Checked all eight combinations of X/Y swapping and axis flips, including coordinate positions, tick direction, axis-label swapping, and preservation of selected cells. Saved each layout and reopened it offline in a fresh browser context; selections, colors, orientation, branch-length mode, and zoom were retained, and further selection worked.
- Checked the English-only controls, removal of the specified headings/notes, the gray `Unselected` legend entry, desktop rendering, and mobile outer-page width. No browser script errors or external resource requests occurred. Adjusted the lower tree margin to keep tip labels above view controls.
- Reviewed both TSV examples, public exports, version metadata, and packaging. Coordinate parsing tests include shuffled rows, numeric-looking IDs, `NA` IDs, duplicate/missing/extra IDs, nonfinite values, and HTML-safe JSON embedding.

Local environment: Python 3.12.14; anndata 0.13.4; Biopython 1.88; NumPy 2.3.5; pandas 3.0.6; SciPy 1.18.1; Matplotlib 3.10.8; pytest 9.1.1. RapidNJ commit: `ed2d36e219d9db16778b941b5054c0fd021b528a`. Python 3.9 and remote CI were not run locally.

## Layout update 0.2.0.4.1

Two local review rounds were completed on 2026-10-20.

- Round 1: reviewed the layout changes and ran the regression suite with real RapidNJ. The original preprocessing, NJ/WLS functions, parameter order, and defaults are unchanged. The existing plotting test now covers all three layouts in both branch-length modes.
- Round 2: installed the built wheel and passed all 32 pytest cases. On the supplied 282-bin dataset, default filtering retained 96 of 535 SNVs. All six layout/branch-length combinations returned the same Newick and WLS diagnostics as the uploaded 0.2.0.4 package.
- Independently compared rectangular plots against the uploaded package. Four combinations of branch-length display and leaf labels/side markers produced identical canvas pixels.
- Checked capitalized layout names, group colors, leaf labels, side markers, and package version alignment. The new layouts use the already declared NumPy and Matplotlib dependencies.

## Previous validation: 0.2.0.4

Two review rounds were completed on 2026-10-20. These are local checks; the included GitHub Actions workflow has not yet run on a published commit.

## Round 1: implementation and regression checks

- 27 pytest cases passed with the real RapidNJ executable.
- Covered dense/sparse inputs, aligned names, zero-depth NaN warnings, rejection of invalid observed ratios/depths, state boundaries, raw-ratio recovery, depth filtering, and exact TSV exclusions.
- Germline decisions matched the benchmark at 2,000 random sites under nine threshold combinations (18,000 decisions).
- E-value threshold decisions were also checked against independent numerical integration.
- An additive eight-cell fixture recovered the known root with WLS powers 0, 1, 2, and 3. Its unrooted splits and pairwise distances were preserved.
- Both plot modes preserved the returned root, topology, and Newick branch lengths; topology-mode tips were aligned.

## Round 2: independent and end-to-end checks

- Independently checked 100 random five-tip WLS cases against bounded numerical optimization, including negative NJ lengths projected to zero. The maximum loss difference was 5.68e-14. Root-to-tip residuals and projected backbone distances were also checked.
- On a 200-cell, 515-site simulated dataset, the package matched the benchmark's retained sites, processed matrix, rooted Newick, and WLS diagnostics in all 16 combinations of germline filtering on/off, two-state/continuous representation, and weight powers 0/1/2/3.
- Checked public API documentation, local documentation links, version metadata, and the supplied six-locus TSV.
- Built source and wheel distributions and verified installation and tests from the wheel.

## Naming and storage update

- Rechecked the renamed API and layer lifecycle: X retains original ratios; assigned states are temporary; tree construction reads the zero-filled output layer.
- 30 pytest cases passed with real RapidNJ, including custom layer names, in-place processing, and rebuilding from the preserved ratios.
- Compared the uploaded package and this update on the same 200-cell, 515-site dataset. Retained sites, processed matrices, rooted Newick, and WLS diagnostics matched in all 32 combinations of germline filtering on/off, two-state/continuous representation, callable depths 1/3, and WLS powers 0/1/2/3.
- Confirmed that germline evidence, RapidNJ execution, WLS optimization, and plotting algorithms are unchanged apart from parameter names.

## Local environment

Python 3.12; anndata 0.10.9; Biopython 1.88; NumPy 2.3.5; SciPy 1.17.0; pandas 2.2.3; Matplotlib 3.10.8; pytest 9.1.1.

RapidNJ commit: `ed2d36e219d9db16778b941b5054c0fd021b528a`.

The upstream repository used as the package starting point was commit `1efc04d95db299b15800230f3e51e8c063e17858` of `RuiyingChenBioinfo/GEI-NJ`.

Python 3.9 is configured in CI but was not locally tested. Version-specific DOI archiving remains to be completed as described in RELEASE.md.
