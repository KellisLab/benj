
from .parse_anndata import setup_args_anndata
from .setup_scanpy import *
from .log import setup_args_log, setup_logging
from .mudata import setup_args_mudata

def parse_args(ap, which=["log"], parse_anndata_prefix=[], load_all:bool=True):
    if "anndata" in which:
        if parse_anndata_prefix is not None and len(parse_anndata_prefix) > 0:
            for pap in parse_anndata_prefix:
                ap = setup_args_anndata(ap, parse_anndata_prefix=pap, load_all=load_all)
        else:
            ap = setup_args_anndata(ap, load_all=load_all)
    if "mudata" in which:
        ap = setup_args_mudata(ap)
    if "scanpy" in which:
        ap = setup_args_scanpy(ap)
    if "log" in which or "logging" in which:
        ap = setup_args_log(ap)
    args = vars(ap.parse_args())
    if "scanpy" in which:
        args = setup_scanpy(**args)
    if "log" in which or "logging" in which:
        args = setup_logging(**args)
    return args
