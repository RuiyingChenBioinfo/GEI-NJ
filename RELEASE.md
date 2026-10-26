# Release GEI-NJ_v0.2.0.5

1. Review the source, tests, documentation, and input encoding. Confirm that `pyproject.toml`, `gei_nj.__version__`, and `CITATION.cff` all report `0.2.0.5`.
2. Install the pinned RapidNJ commit shown in the README. Run:

   ```bash
   python -m pip install -e ".[test]" build
   python -m pytest
   python -m build
   ```

   Confirm that both distributions include the coordinate-viewer HTML template. Generate a viewer with `view_tree_coords()`, check node selection, Hide tip names, and coordinate flips/swapping in a browser, and reopen an exported HTML to check that its selection, colors, and tip-name visibility persist.

3. Commit the reviewed files and run the GitHub Actions workflow on that commit. Confirm both Python jobs pass. Record the commit used for the paper.
4. Create the `GEI-NJ_v0.2.0.5` tag at that exact commit and publish its GitHub release. Attach the wheel and source distribution built from that commit.
5. Archive that exact tagged source in Zenodo, either through its GitHub integration or a manual software deposit. Include the GitHub tag, commit, software version, and RapidNJ commit in the record. Cite the **version-specific DOI** in the paper.

If the DOI must be embedded in the release itself, reserve it in a Zenodo draft first. Add the reserved DOI to `CITATION.cff` and the README before the final commit, rerun CI, and archive that commit. Do not move a published tag to add metadata.

The workflow uses RapidNJ commit `ed2d36e219d9db16778b941b5054c0fd021b528a`. Python dependencies are declared in `pyproject.toml`; record the resolved environment alongside the paper's analysis outputs.
