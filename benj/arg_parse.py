
from .parse_anndata import setup_args_anndata
from .setup_scanpy import *
from .log import setup_args_log, setup_logging

def parse_args(ap, which=["log"]):
    if "anndata" in which:
        ap = setup_args_anndata(ap)
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
