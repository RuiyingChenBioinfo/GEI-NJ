def comb_depth_ratio(ad_depth, ad_ratio):
    """
    Combine SNV depth (ad_depth) and ratio (ad_ratio) into a single AnnData.

    Parameters
    ----------
    ad_depth : AnnData
        Depth data. Must cover the same obs/var as ad_ratio.
    ad_ratio : AnnData
        SNV ratio data. Values must be in [0, 1].

    Returns
    -------
    AnnData
        Copy of ad_ratio with depth stored in layers["depth"].
    """
    import numpy as np
    from scipy import sparse

    X = ad_ratio.X
    if sparse.issparse(X):
        data = X.data
    else:
        data = X

    if np.nanmin(data) < 0 or np.nanmax(data) > 1:
        raise ValueError("Ratio data should be in [0,1]")

    adata = ad_ratio.copy()
    r = ad_depth[adata.obs_names, adata.var_names]
    adata.layers["depth"] = r.X

    return adata


def filt_snv_by_depth_ratio(
    adata,
    min_tot_depth=20,
    min_mut_depth=0,
    min_mut_avg_ratio=0,
    x_layer=None,
    depth_layer="depth",
):
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    # checks before running
    if depth_layer not in adata.layers:
        raise KeyError(f"Missing '{depth_layer}' in adata.layers.")
    depth = adata.layers[depth_layer]

    X = adata.layers[x_layer] if x_layer is not None else adata.X

    # depth could be per-cell (n_obs,) or (n_obs, 1), or per-cell-per-var (n_obs, n_vars)
    per_cell = (getattr(depth, "ndim", 2) == 1) or (hasattr(depth, "shape") and len(depth.shape) == 2 and depth.shape[1] == 1)

    # total depth per SNV
    if per_cell:
        d = np.asarray(depth.toarray()).ravel() if sp.issparse(depth) else np.asarray(depth).ravel()
        denom_scalar = float(d.sum())
        # treat D_ij = d_i for all j (consistent with weighted_per_var per_cell branch)
        snv_depth_sums = pd.Series(np.full(adata.n_vars, denom_scalar, dtype=float), index=adata.var_names)
    else:
        if sp.issparse(depth):
            snv_depth_arr = np.asarray(depth.sum(axis=0)).ravel()
        else:
            snv_depth_arr = np.asarray(depth).sum(axis=0)

        if snv_depth_arr.shape[0] != adata.n_vars:
            raise ValueError("Depth matrix's number of columns does not match adata.var_names.")

        snv_depth_sums = pd.Series(snv_depth_arr.astype(float, copy=False), index=adata.var_names)

    # mutated reads sums: sum_i X_ij * D_ij  (only compute if needed for min_mut_depth or min_mut_avg_ratio)
    need_num = (float(min_mut_depth) > 0) or (float(min_mut_avg_ratio) > 0)

    if need_num:
        if per_cell:
            d = np.asarray(depth.toarray()).ravel() if sp.issparse(depth) else np.asarray(depth).ravel()
            if sp.issparse(X):
                num = np.asarray(X.multiply(d[:, None]).sum(axis=0)).ravel()
            else:
                num = (np.asarray(X) * d[:, None]).sum(axis=0)
            num = num.astype(float, copy=False)
        else:
            if sp.issparse(X) and sp.issparse(depth):
                num = np.asarray(X.multiply(depth).sum(axis=0)).ravel()
            elif sp.issparse(X) and not sp.issparse(depth):
                # sparse X elementwise-multiply dense depth is supported
                num = np.asarray(X.multiply(np.asarray(depth)).sum(axis=0)).ravel()
            elif (not sp.issparse(X)) and sp.issparse(depth):
                num = np.asarray(depth.multiply(np.asarray(X)).sum(axis=0)).ravel()
            else:
                num = (np.asarray(X) * np.asarray(depth)).sum(axis=0)
            num = num.astype(float, copy=False)

        if num.shape[0] != adata.n_vars:
            raise ValueError("Mutated reads matrix's number of columns does not match adata.var_names.")

        snv_mut_sums = pd.Series(num, index=adata.var_names)
    else:
        snv_mut_sums = pd.Series(np.zeros(adata.n_vars, dtype=float), index=adata.var_names)

    # weighted average ratio per SNV (same as weighted_per_var)
    if float(min_mut_avg_ratio) > 0:
        if per_cell:
            denom = float(snv_depth_sums.iloc[0])  # scalar
            mut_avg_ratio = np.divide(
                snv_mut_sums.to_numpy(),
                denom,
                out=np.full(adata.n_vars, np.nan, dtype=float),
                where=(denom != 0),
            )
        else:
            denom = snv_depth_sums.to_numpy()
            mut_avg_ratio = np.divide(
                snv_mut_sums.to_numpy(),
                denom,
                out=np.full(adata.n_vars, np.nan, dtype=float),
                where=(denom != 0),
            )
        mut_avg_ratio = pd.Series(mut_avg_ratio, index=adata.var_names)
    else:
        # not used in mask
        mut_avg_ratio = None

    # build mask
    mask = (snv_depth_sums >= float(min_tot_depth))

    if float(min_mut_depth) > 0:
        mask &= (snv_mut_sums >= float(min_mut_depth))

    if float(min_mut_avg_ratio) > 0:
        mask &= (mut_avg_ratio >= float(min_mut_avg_ratio))

    snv_to_keep = snv_depth_sums.index[mask]

    print(f"SNV number before filter: {adata.n_vars}")
    print(f"SNV number after filter: {len(snv_to_keep)}")

    return adata[:, snv_to_keep].copy()


def filt_adata_by_mutation_nnz(
    adata,
    min_nnz: int = 2,
    binarize_nonzero_to_one: bool = False,
    to_csr: bool = True,
    verbose: bool = True,
):
    """
    Filter AnnData columns (mutations) by the number of non-zero entries per column,
    and optionally binarize all non-zero values to 1.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object (cells x mutations) stored in adata.X.
    min_nnz : int, default=2
        Keep columns that have at least `min_nnz` non-zero entries.
    binarize_nonzero_to_one : bool, default=False
        If True, convert any non-zero value (e.g., 0.5 or 1) to 1.
        Note: for sparse matrices, this sets ALL stored data entries to 1.0,
        matching the original behavior (including the rare case of explicitly stored zeros).
    to_csr : bool, default=True
        If True and the resulting matrix is sparse, convert it to CSR format.
    verbose : bool, default=True
        If True, print basic stats and how many columns are kept.

    Returns
    -------
    adata_filt : AnnData
        A filtered copy of the input AnnData.
    """
    import numpy as np
    import scipy.sparse as sp

    X = adata.X
    n_cells, n_muts = adata.n_obs, adata.n_vars

    # 1) Compute per-column nnz on the original matrix (no full copy)
    if sp.issparse(X):
        # Fast paths for common formats
        if sp.isspmatrix_csr(X):
            # Count occurrences of each column index in CSR indices
            nnz_per_gene = np.bincount(X.indices, minlength=n_muts)
        elif sp.isspmatrix_csc(X):
            # In CSC, indptr encodes nnz per column
            nnz_per_gene = np.diff(X.indptr)
        else:
            nnz_per_gene = np.asarray(X.getnnz(axis=0)).ravel()

        nonzero_total = int(X.nnz)
    else:
        nnz_per_gene = np.count_nonzero(X, axis=0)
        nonzero_total = int(np.count_nonzero(X))

    keep_mask = nnz_per_gene >= int(min_nnz)
    n_keep = int(keep_mask.sum())

    if verbose:
        total_entries = n_cells * n_muts
        ratio = nonzero_total / total_entries if total_entries else float("nan")
        print(f"Input shape: {n_cells} cells x {n_muts} mutations")
        print(f"Non-zero ratio: {ratio:.6g}")
        print(f"Number of kept mutations (nnz>={min_nnz}): {n_keep}")

        # Same optional diagnostics as your original implementation
        if sp.issparse(X):
            vals = X.data
            count_1 = int(np.count_nonzero(np.isclose(vals, 1.0)))
            count_05 = int(np.count_nonzero(np.isclose(vals, 0.5)))
            print(f"Original sparse data counts: 1.0 -> {count_1}, 0.5 -> {count_05}")

    # 2) Slice AnnData once (avoid double slicing X)
    adata_filt = adata[:, keep_mask].copy()

    # 3) Optional binarization on the filtered submatrix only
    if binarize_nonzero_to_one:
        if sp.issparse(adata_filt.X):
            # Match original behavior: set ALL stored data entries to 1.0
            adata_filt.X.data[:] = 1.0
        else:
            adata_filt.X = (adata_filt.X != 0).astype(np.float32)

    # 4) Convert to CSR only if requested and needed
    if to_csr and sp.issparse(adata_filt.X) and (not sp.isspmatrix_csr(adata_filt.X)):
        adata_filt.X = adata_filt.X.tocsr()

    return adata_filt


def revise_snv_ratio(
    adata,
    layer_name="original_ratio",
    low=0.1,
    high=0.9,
    mid=0.5,
    copy=True,
    dtype="float32",
):
    """
    1) Copy adata.X into a new layer (default: "original_ratio").
    2) Discretize adata.X in-place (or on a copy) using configurable thresholds:
         - values < low   -> 0
         - values > high  -> 1
         - otherwise      -> mid
       Note: values exactly equal to low/high fall into the "otherwise" bucket.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    layer_name : str
        Name of the layer that will store the original (pre-discretization) X.
    low : float
        Lower threshold.
    high : float
        Upper threshold.
    mid : float
        Value assigned to entries within [low, high] (inclusive).
    copy : bool
        If True, operate on a copy and return it; if False, modify the input adata in-place.
    dtype : str or numpy dtype
        Output dtype for dense X and the stored original layer (dense). Sparse keeps its dtype.

    Returns
    -------
    AnnData
        AnnData containing layers[layer_name] = original X, and modified X.
    """
    import numpy as np
    from scipy import sparse

    if low > high:
        raise ValueError(f"`low` must be <= `high`. Got low={low}, high={high}.")

    ad = adata.copy() if copy else adata

    X = ad.X

    # Store the original X into a new layer
    if sparse.issparse(X):
        ad.layers[layer_name] = X.copy()
        M = X.tocsr(copy=True)
        data = M.data

        mask_low = data < low
        mask_high = data > high
        mask_mid = (~mask_low) & (~mask_high)

        data[mask_low] = 0.0
        data[mask_high] = 1.0
        data[mask_mid] = float(mid)

        M.eliminate_zeros()
        ad.X = M
    else:
        ad.layers[layer_name] = np.array(X, copy=True).astype(dtype, copy=False)

        arr = np.array(X, copy=True)
        out = np.full(arr.shape, float(mid), dtype=dtype)
        out[arr < low] = 0.0
        out[arr > high] = 1.0

        ad.X = out

    return ad


def add_callable_layer(
    adata,
    min_callable_depth=3,
    depth_layer="depth",
    callable_layer="callable",
    copy=True,
    dtype="int8",
    verbose=True,
):
    """
    Add a binary callable layer to an AnnData object.

    A CMB x SNV entry is considered callable if its sequencing depth is greater
    than or equal to min_callable_depth.
    """
    import numpy as np
    import scipy.sparse as sp

    if depth_layer not in adata.layers:
        raise KeyError(f"Missing '{depth_layer}' in adata.layers.")

    ad = adata.copy() if copy else adata
    depth = ad.layers[depth_layer]

    if sp.issparse(depth):
        callable_mat = depth.copy().tocsr()
        callable_mat.data = (callable_mat.data >= min_callable_depth).astype(dtype)
        callable_mat.eliminate_zeros()

        callable_total = int(callable_mat.sum())
        n_total = callable_mat.shape[0] * callable_mat.shape[1]

    else:
        depth_arr = np.asarray(depth)

        if depth_arr.ndim == 1:
            depth_arr = depth_arr[:, None]

        if depth_arr.shape[1] == 1 and ad.n_vars > 1:
            callable_mat = np.repeat(
                depth_arr >= min_callable_depth,
                ad.n_vars,
                axis=1,
            ).astype(dtype)
        else:
            if depth_arr.shape != ad.shape:
                raise ValueError(
                    f"Depth layer shape {depth_arr.shape} does not match adata shape {ad.shape}."
                )

            callable_mat = (depth_arr >= min_callable_depth).astype(dtype)

        callable_total = int(callable_mat.sum())
        n_total = callable_mat.shape[0] * callable_mat.shape[1]

    ad.layers[callable_layer] = callable_mat

    global_callable_fraction = callable_total / n_total if n_total > 0 else float("nan")

    if verbose:
        print(f"Callable layer added to adata.layers['{callable_layer}']")
        print(f"Minimum callable depth: {min_callable_depth}")
        print(f"Total SNV entries: {n_total}")
        print(f"Callable SNV entries: {callable_total}")
        print(f"Global callable fraction: {global_callable_fraction:.4f}")

    return ad


def impute_non_callable_by_snv_mean(
    adata,
    callable_layer="callable",
    output_layer="geinj_state_imputed",
    copy=True,
    dtype="float32",
    verbose=True,
):
    """
    Impute non-callable entries in revised SNV state matrix by the per-SNV mean
    calculated from callable entries.

    This function assumes adata.X contains revised SNV states, such as 0, 0.5,
    and 1, and adata.layers[callable_layer] contains binary callability.
    """
    import numpy as np
    import scipy.sparse as sp

    if callable_layer not in adata.layers:
        raise KeyError(f"Missing adata.layers['{callable_layer}'].")

    ad = adata.copy() if copy else adata

    X = ad.X.toarray() if sp.issparse(ad.X) else np.asarray(ad.X)
    C = ad.layers[callable_layer]
    C = C.toarray() if sp.issparse(C) else np.asarray(C)

    X = X.astype(dtype, copy=True)
    C = C.astype(bool, copy=False)

    if X.shape != C.shape:
        raise ValueError(
            f"Shape mismatch: adata.X has shape {X.shape}, "
            f"but callable layer has shape {C.shape}."
        )

    X_masked = X.copy()
    X_masked[~C] = np.nan

    snv_mean = np.nanmean(X_masked, axis=0)

    global_mean = np.nanmean(X_masked)
    snv_mean = np.where(np.isnan(snv_mean), global_mean, snv_mean)

    missing_i, missing_j = np.where(~C)
    X_masked[missing_i, missing_j] = snv_mean[missing_j]

    X_imputed = X_masked.astype(dtype, copy=False)

    ad.layers[output_layer] = X_imputed
    ad.X = X_imputed

    if verbose:
        n_total = X.size
        n_non_callable = int((~C).sum())
        n_callable = int(C.sum())
        print(f"Non-callable entries imputed by per-SNV callable mean")
        print(f"Callable entries: {n_callable}")
        print(f"Non-callable entries: {n_non_callable}")
        print(f"Total SNV entries: {n_total}")
        print(f"Imputed matrix stored in adata.layers['{output_layer}']")
        print(f"adata.X was updated to the imputed matrix")

    return ad
