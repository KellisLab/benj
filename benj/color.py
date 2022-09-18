

def cm_benj(base="viridis", count=100, zero=[0.4, 0.4, 0.4, 1]):
    from matplotlib.cm import get_cmap
    import matplotlib.colors as mcolors
    import numpy as np
    cm = get_cmap(base)(np.linspace(0, 1, count-1))
    colors = np.vstack((zero, cm))
    return mcolors.LinearSegmentedColormap.from_list('benj', colors)
