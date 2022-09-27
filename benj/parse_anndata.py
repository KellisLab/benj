
def setup_args_anndata(ap):
    ap.add_argument("-a", "--annotation", default="")
    ap.add_argument("--min", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--max", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subset", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subsample", default=1, type=int)
    ap.add_argument("--split", default=",")
    return ap

def parse_anndata(h5ad, **args):
    import os
    import numpy as np
    if os.path.exists(h5ad):
        import anndata
        adata = anndata.read(args["input"], backed="r")
    else:
        import sys
        raise RuntimeError("File \"%s\" does not exist." % args["input"])
    obs = adata.obs.copy()
    flag = np.repeat(True, obs.shape[0])
    if args.get("annotation") and os.path.exists(args["annotation"]):
        import pandas as pd
        annot = pd.read_csv(args["annotation"], index_col=0, sep="\t")
        obs = obs.loc[annot.index.values,:].copy()
        for cn in annot.columns.values:
            obs[cn] = annot[cn].values
    if args.get("subset") and isinstance(args["subset"], list):
        for ss in args["subset"]:
            S = [s.strip() for s in ss.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].isin(S[1].split(args["split"])).values
                if np.sum(flag) == 0:
                    raise RuntimeError("%s not in %s" % (S[1].split(args["split"]), S[0]))
                else:
                    print(S[0], "is not a valid column")
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
    return adata[obs.index.values[flag], :].to_memory()
