
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

def setup_args_anndata(ap):
    ap.add_argument("-a", "--annotation", nargs="+")
    ap.add_argument("--min", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--max", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subset", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subsample", default=1, type=int)
    ap.add_argument("--split", default=",")
    return ap

def parse_anndata(**args):
    import os
    import logging
    import numpy as np
    h5ad = args.get("h5ad", args.get("input", ""))
    if os.path.exists(h5ad):
        import anndata
        adata = anndata.read(h5ad, backed="r")
    else:
        import sys
        raise RuntimeError("File \"%s\" does not exist." % h5ad)
    obs = adata.obs.copy()
    flag = np.repeat(True, obs.shape[0])
    if args.get("annotation"):
        import pandas as pd
        for fname in args["annotation"]:
            if os.path.exists(fname):
                annot = pd.read_csv(fname, index_col=0, sep="\t", low_memory=False)
                annot = annot.loc[annot.index.isin(obs.index.values), :]
                obs = obs.loc[annot.index.values,:].copy()
                for cn in annot.columns.values:
                    obs[cn] = annot[cn].values
            else:
                print("Warning: File", fname, "does not exist")
        flag = np.repeat(True, obs.shape[0])
    if args.get("subset") and isinstance(args["subset"], list):
        for ss in args["subset"]:
            S = [s.strip() for s in ss.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                from ast import literal_eval
                sval = S[1].split(args["split"])
                try:
                    sval = [literal_eval(x) for x in sval]
                except ValueError:
                    pass
                dt = obs[S[0]].dtype
                if dt.name == "category":
                    dt = str
                sval = np.asarray(sval, dtype=dt)
                flag &= obs[S[0]].isin(sval).values
                if np.sum(flag) == 0:
                    raise RuntimeError("%s not in %s" % (S[1].split(args["split"]), S[0]))
    if args.get("min") and isinstance(args["min"], list):
        for kv in args["min"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].values >= float(S[1])
            else:
                print(S[0], "is not a valid column")
    if args.get("max") and isinstance(args["max"], list):
        for kv in args["max"]:
            S = [s.strip() for s in kv.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].values <= float(S[1])
            else:
                print(S[0], "is not a valid column")
    if args.get("subsample") and args["subsample"] > 1:
        ss = args["subsample"]
        flag = np.ravel(np.where(flag))
        flag = np.random.choice(flag, len(flag)//args["subsample"], replace=False)
    adata = adata[obs.index.values[flag], :]
    if args.get("backed") is None:
        adata = adata.to_memory()
    for cn in np.setdiff1d(obs.columns, adata.obs.columns):
        adata.obs[cn] = obs[cn]
    logger = logging.getLogger(args.get("logger", "benj"))
    logging.basicConfig(level=logging.INFO)
    logger.info("Read %d cells" % adata.shape[0])
    return adata
