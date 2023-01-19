
def run_nmf(X, n_components, alpha=1.0, max_iter=200, verbose=True, random_state=42, init="nndsvd"):
    from scopen.MF import NMF
    model = NMF(n_components=n_components,
                random_state=random_state,
                init=init,
                alpha=alpha,
                l1_ratio=0,
                max_iter=max_iter,
                verbose=verbose)
    w_hat = model.fit_transform(X=X)
    h_hat = model.components_
    print(f"ranks: {n_components}, fitting error: {model.reconstruction_err_}")
    return [w_hat, h_hat, model.reconstruction_err_]

def _run_nmf(arguments):
    X, n_components, alpha, max_iter, verbose, random_state, init = arguments
    return run_nmf(X=X, n_components=n_components, alpha=alpha, max_iter=max_iter, verbose=verbose, random_state=random_state, init=init)


def nmf(X, a_min=2, a_max=30, a_step=1, alpha=1.0, init="nndsvd", random_state=42, verbose=True, max_iter=200, nproc=16):
    import numpy as np
    from multiprocessing import Pool
    from kneed import KneeLocator
    from functools import partial
    arguments_list = list()
    n_components_list = np.arange(a_min, a_max+1, a_step)
    w_hat_dict = {}
    h_hat_dict = {}
    error_list = []
    for n_components in n_components_list:
        arguments = (X, n_components, alpha, max_iter, verbose, random_state, init)
        arguments_list.append(arguments)
    with Pool(processes=nproc) as pool:
        res = pool.map(run_nmf, arguments_list)
    print("res:")
    print(res)
    for i, n_components in enumerate(n_components_list):
            w_hat_dict[n_components] = res[i][0]
            h_hat_dict[n_components] = res[i][1]
            error_list.append(res[i][2])
    if len(error_list) > 1:
        kl = KneeLocator(n_components_list, error_list,
                         S=1.0, curve="convex", direction="decreasing")
        return w_hat_dict[kl.knee], h_hat_dict[kl.knee]
    else:
        return w_hat_dict[n_components_list[0]], h_hat_dict[n_components_list[0]]
