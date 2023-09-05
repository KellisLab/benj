#!/usr/bin/env python3

def run(**kwargs):
    from benj.count_atac import count_atac
    count_atac(**kwargs)
    
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fragments", required=True)
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("--cell-metadata", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--peaks-bed", dest="bed", action="store_true")
    ap.add_argument("--peaks-tsv", dest="bed", action="store_false")
    ap.add_argument("-b", "--blacklist", default=None)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--max-value", type=int, default=127)
    ap.add_argument("--stranded", dest="stranded", action="store_true")
    ap.add_argument("--unstranded", dest="stranded", action="store_false")
    ap.add_argument("--qc", nargs="+", default=["peakType"])
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(bed=False, stranded=False)
    args = vars(ap.parse_args())
    run(**args)
