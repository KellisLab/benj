

def window(var, interval="interval", size:int=10000000, step:int=500000):
        """Returns a matrix transformation, var x N shape.

        Must use X * W[,i]
        import p"""
        import numpy as np
        import pandas as pd
        import scipy.sparse
        import pyranges
        from .utils import as_ranges
        gr = as_ranges(var, interval=interval)
        bounds = gr.apply(lambda x: (x.Start.min(), x.End.max()), as_pyranges=False)
        wf = pd.DataFrame(bounds, index=["Start", "End"]).T.reset_index(names="Chromosome")
        out = []
        for i in range(0, size, step):
                out.append(wf.copy())
                wf["Start"] += step
        wf = pyranges.PyRanges(pd.concat(out)).window(size)
        wf = wf.sort().df.drop_duplicates()
        wf["_index"] = np.arange(len(wf))
        wf = pyranges.PyRanges(wf)
        pf = gr.join(wf).df.loc[:, ["_index", "_index_b"]]
        row = pd.Index(var[interval]).get_indexer(pf["_index"].values)
        col = pf["_index_b"].values
        S = scipy.sparse.csr_matrix((np.ones(len(row), dtype="i1"), (row, col)),
                                    shape=(var.shape[0], len(wf)), dtype="i1")
        return S


def window_transform_full(X, W):
        import numpy as np
        import scipy.sparse
        R = np.ones((X.shape[0], 1), dtype="i1")
        R = scipy.sparse.csr_matrix(R)
        return scipy.sparse.kron(R, X).dot(W)
