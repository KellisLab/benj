
def setup_args_mudata(ap):
    ap.add_argument("-a", "--annotation", nargs="+")
    ap.add_argument("--min", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--max", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subset", nargs="+", metavar="KEY=VALUE")
    ap.add_argument("--subsample", default=1, type=int)
    ap.add_argument("--split", default=",")
    ap.add_argument("--intersect-obs", action="store_true", dest="intersect_obs")
    ap.add_argument("--no-intersect-obs", action="store_false", dest="intersect_obs")
    ap.set_defaults(intersect_obs=False)
    return ap

def parse_mudata(**args):
    import os
    import numpy as np
    h5mu = args.get("h5mu", args.get("input", ""))
    if os.path.exists(h5mu):
        import muon
        mdata = muon.read_h5mu(h5mu, backed="r")
    else:
        import sys
        raise RuntimeError("File \"%s\" does not exist." % h5mu)
    obs = mdata.obs.copy()
    flag = np.repeat(True, obs.shape[0])
    if args.get("annotation"):
        import pandas as pd
        for fname in args["annotation"]:
            if os.path.exists(fname):
                annot = pd.read_csv(fname, index_col=0, sep="\t")
                obs = obs.loc[annot.index.values,:].copy()
                for cn in annot.columns.values:
                    obs[cn] = annot[cn].values
            else:
                print("Warning: File", fname, "does not exist")
    if args.get("subset") and isinstance(args["subset"], list):
        for ss in args["subset"]:
            S = [s.strip() for s in ss.split("=")]
            if len(S) == 2 and S[0] in obs.columns:
                flag &= obs[S[0]].isin(S[1].split(args["split"])).values
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
    mdata = mdata[obs.index.values[flag], :].to_memory()
    for cn in np.setdiff1d(obs.columns, mdata.obs.columns):
        mdata.obs[cn] = obs[cn]
    if mdata.filename is not None:
        mdata.file._to_memory_mode()
    if args.get("intersect_obs"):
        import muon
        muon.pp.intersect_obs(mdata)
    return mdata
