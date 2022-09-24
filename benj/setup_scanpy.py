
from .color import cm_benj

def setup_args_scanpy(ap):
    ap.add_argument("-f", "--figdir", default="figures")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--frameon", dest="frameon", action="store_true")
    ap.add_argument("--no-frameon", dest="frameon", action="store_false")
    ap.add_argument("--verbosity", type=int, default=2)
    ap.set_defaults(frameon=True)
    return ap

def setup_scanpy(**args):
    import scanpy as sc
    if args.get("figdir"):
        sc.settings.figdir = args["figdir"]
    cm = cm_benj()
    if args.get("dpi"):
        dpi = args["dpi"]
    else:
        dpi = 600
    frameon = True
    if args.get("frameon") and isinstance(args["frameon"], bool):
        frameon = args["frameon"]
    if args.get("verbosity") and isinstance(args["verbosity"], int):
        sc.settings.verbosity = args["verbosity"]
    else:
        sc.settings.verbosity = 2
    sc.set_figure_params(dpi_save=dpi, color_map=cm_benj(), frameon=frameon)
    return args
