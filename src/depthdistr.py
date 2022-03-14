import matplotlib.pyplot as plt
import numpy as np

def sample_depths(nevents, mean_depth = 20, isaftershock=False, doplot=False):
    # depths are assumed to be lognormally distrbuted...
    mu = np.log(mean_depth)
    if isaftershock:
        sigma = 0.1
    else:
        sigma = 0.2
    d = np.random.lognormal(mean=mu, sigma=sigma, size=nevents)
    if doplot:
        fig, ax = plt.subplots(1, 1)
        ax.hist(d, density=True)
        ax.set_xlabel('depth (km)')
    return d

