
def convert_dict(dct):
    import anndata
    from collections import OrderedDict
    if type(dct) not in [OrderedDict, anndata.compat._overloaded_dict.OverloadedDict]:
        return dct
    return {k: convert_dict(v) for k, v in dct.items()}

def remove_ordereddict(mdata):
    for k in mdata.mod.keys():
        mdata.mod[k].uns = convert_dict(mdata.mod[k].uns)
    return mdata

def mu_load_with_aggr(feature_matrix_h5, aggr_csv, library_id="library_id", *extra_cols):
    from .timer import template as stopwatch
    import pandas as pd
    import muon as mu
    af = pd.read_csv(aggr_csv)
    if sw is None:
        sw = stopwatch()
    with sw("Loading H5") as _:
        mdata = mu.read_10x_h5(h5)
    mdata.var_names_make_unique()
    mdata = remove_ordereddict(mdata)
    ### add aggr data
    idx = mdata.obs_names.str.replace(r"^[^-]+-", "", regex=True).astype(int).values - 1
    mdata.obs[library_id] = af["library_id"].values[idx]
    for col in extra_cols:
        mdata.obs[col] = af[col].values[idx]
    for k in mdata.mod.keys():
        with sw("Converting .X to int for .mod[%s]" % k) as _:
            mdata.mod[k].X = mdata.mod[k].X.astype(int)
        for col in list(extra_cols) + [library_id]:
            mdata.mod[k].obs[col] = mdata.obs[col]
    return mdata
