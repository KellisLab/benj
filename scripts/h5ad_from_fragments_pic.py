#!/usr/bin/env python3

def run(**kwargs):
    from benj import pic_create_h5ad
    pic_create_h5ad(**kwargs)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fragments", required=True)
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("-o", "--output", required=True, dest="h5ad")
    ap.add_argument("-p", "--peaks", required=True, dest="peak_bed")
    ap.add_argument("--tss", default=None, dest="tss_bed")
    ap.add_argument("-b", "--blacklist", default=None, dest="blacklist_bed")
    ap.add_argument("--compression", type=int, default=9)
    ap.add_argument("--sample-column", default="Sample", dest="sample_name")
    ap.add_argument("--extend", type=int, default=5)
    ap.add_argument("--promoter-upstream", type=int, default=2000)
    ap.add_argument("--promoter-downstream", type=int, default=100)
    args = vars(ap.parse_args())
    run(**args)
