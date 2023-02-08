def remove_ordereddict(mdata):
    def convert_dict(dct):
        import anndata
        from collections import OrderedDict
        if type(dct) not in [OrderedDict, anndata.compat._overloaded_dict.OverloadedDict]:
            return dct
        return {k: convert_dict(v) for k, v in dct.items()}
    for k in mdata.mod.keys():
        mdata.mod[k].uns = convert_dict(mdata.mod[k].uns)
    return mdata
