
def annotation_to_numbered(df, sample_order, index_unique="#", startnum:int=1, label=None):
    import pandas as pd
    out = []
    for i, sample in enumerate(sample_order):
        flag = df.index.str.contains("^[ACGT]+-1%s%s$" % (index_unique, sample))
        nf = df.loc[flag, :].copy()
        nf.index = nf.index.str.replace("1%s%s$" % (index_unique, sample),
                                        "%d" % (startnum + i),
                                        regex=True)
        if label is not None:
            nf[label] = sample
        out.append(nf)
    return pd.concat(out)

def setup_args_anndata(ap, parse_anndata_prefix="", load_all:bool=True):
    load_def = []
    if load_all:
        load_def = ["all"]
    if parse_anndata_prefix:
        parse_anndata_prefix = parse_anndata_prefix.replace("_", "-") + "-"
        ap.add_argument("--%sannotation" % parse_anndata_prefix, nargs="+")
        ap.add_argument("--%smin" % parse_anndata_prefix, nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--%smax" % parse_anndata_prefix, nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--%ssubset" % parse_anndata_prefix, nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--%ssubsample" % parse_anndata_prefix, default=1, type=int)
        ap.add_argument("--%ssplit" % parse_anndata_prefix, default=",")
        ap.add_argument("--%sload" % parse_anndata_prefix, nargs="+", default=load_def, help="Load X, layers, or all, or none")
    else:
        ap.add_argument("-a", "--annotation", nargs="+")
        ap.add_argument("--min", nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--max", nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--subset", nargs="+", metavar="KEY=VALUE")
        ap.add_argument("--subsample", default=1, type=int)
        ap.add_argument("--split", default=",")
        ap.add_argument("--load", nargs="+", default=load_def, help="Load X, layers, or all, or none")
    return ap

def parse_anndata_options(parse_anndata_prefix="", **args):
    """ Use argparse derived args and split into arguments to aggregate_concat"""
    if parse_anndata_prefix:
        parse_anndata_prefix = parse_anndata_prefix.replace("-", "_") + "_"
    h5ad = args.get(parse_anndata_prefix + "h5ad",
                    args.get(parse_anndata_prefix + "input", ""))
    ### h5ad can be either str or list[str]
    load = args.get(parse_anndata_prefix + "load", "all")
    out = {"h5ad": h5ad, "load": load}
    if args.get(parse_anndata_prefix + "subset") and isinstance(args[parse_anndata_prefix + "subset"], list):
        subset = {}
        for ss in args[parse_anndata_prefix + "subset"]:
            S = [s.strip() for s in ss.split("=")]
            if len(S) == 2:
                from ast import literal_eval
                sval = S[1].split(args.get(parse_anndata_prefix + "split", ","))
                try:
                    sval = [literal_eval(x) for x in sval]
                except:
                    pass
                subset[S[0]] = sval
        out["subset"] = subset
    if args.get(parse_anndata_prefix + "min") and isinstance(args[parse_anndata_prefix + "min"], list):
        obs_min = {}
        for kv in args[parse_anndata_prefix + "min"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2:
                obs_min[S[0]] = float(S[1])
        out["obs_min"] = obs_min
    if args.get(parse_anndata_prefix + "max") and isinstance(args[parse_anndata_prefix + "max"], list):
        obs_max = {}
        for kv in args[parse_anndata_prefix + "max"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2:
                obs_max[S[0]] = float(S[1])
        out["obs_max"] = obs_max
    if args.get(parse_anndata_prefix + "annotation"):
        out["obs_annotation"] = args["annotation"]
    return out

def parse_anndata(parse_anndata_prefix="", **args):
    """Legacy option, uses aggregate_concat and aggregate_load"""
    from .aggregate import aggregate_concat, aggregate_load
    opt = parse_anndata_options(parse_anndata_prefix=parse_anndata_prefix, **args)
    adata = aggregate_concat(**opt)
    if opt.get("load"):
        return aggregate_load(adata, which=opt["load"])
    else:
        return adata

def _old_parse_anndata(parse_anndata_prefix="", **args):
    import os
    import logging
    import numpy as np
    if parse_anndata_prefix:
        parse_anndata_prefix = parse_anndata_prefix.replace("-", "_") + "_"
    h5ad = args.get(parse_anndata_prefix + "h5ad",
                    args.get(parse_anndata_prefix + "input", ""))
    if os.path.exists(h5ad):
        import anndata
        adata = anndata.read(h5ad, backed="r")
    else:
        import sys
        raise RuntimeError("File \"%s\" does not exist." % h5ad)
    obs = adata.obs.copy()
    flag = np.repeat(True, obs.shape[0])
    if args.get(parse_anndata_prefix + "annotation"):
        import pandas as pd
        for fname in args[parse_anndata_prefix + "annotation"]:
            if os.path.exists(fname):
                annot = pd.read_csv(fname, index_col=0, sep="\t", low_memory=False)
                annot = annot.loc[annot.index.isin(obs.index.values), :]
                obs = obs.loc[annot.index.values,:].copy()
                for cn in annot.columns.values:
                    obs[cn] = annot[cn].values
            else:
                print("Warning: File", fname, "does not exist")
        flag = np.repeat(True, obs.shape[0])
    if args.get(parse_anndata_prefix + "subset") and isinstance(args[parse_anndata_prefix + "subset"], list):
        for ss in args[parse_anndata_prefix + "subset"]:
            S = [s.strip() for s in ss.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                from ast import literal_eval
                sval = S[1].split(args[parse_anndata_prefix + "split"])
                try:
                    sval = [literal_eval(x) for x in sval]
                except:
                    pass
                dt = obs[S[0]].dtype
                if dt.name == "category":
                    dt = str
                sval = np.asarray(sval, dtype=dt)
                flag &= obs[S[0]].isin(sval).values
                if np.sum(flag) == 0:
                    raise RuntimeError("%s not in %s" % (S[1].split(args[parse_anndata_prefix + "split"]), S[0]))
    if args.get(parse_anndata_prefix + "min") and isinstance(args[parse_anndata_prefix + "min"], list):
        for kv in args[parse_anndata_prefix + "min"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].values >= float(S[1])
            else:
                print(S[0], "is not a valid column")
    if args.get(parse_anndata_prefix + "max") and isinstance(args[parse_anndata_prefix + "max"], list):
        for kv in args[parse_anndata_prefix + "max"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].values <= float(S[1])
            else:
                print(S[0], "is not a valid column")
    if args.get(parse_anndata_prefix + "subsample") and args[parse_anndata_prefix + "subsample"] > 1:
        ss = args[parse_anndata_prefix + "subsample"]
        flag = np.ravel(np.where(flag))
        flag = np.random.choice(flag, len(flag)//args[parse_anndata_prefix + "subsample"], replace=False)
    adata = adata[obs.index.values[flag], :]
    if args.get(parse_anndata_prefix + "backed") is None or not args.get(parse_anndata_prefix + "backed"):
        adata = adata.to_memory()
    for cn in np.setdiff1d(obs.columns, adata.obs.columns):
        adata.obs[cn] = obs[cn]
    logger = logging.getLogger(args.get("logger", "benj"))
    logging.basicConfig(level=logging.INFO)
    logger.info("Read %d cells" % adata.shape[0])
    return adata
