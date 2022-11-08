import matplotlib.pyplot as plt
import numpy as np
from random import choices, seed

def distr_empirical(M, minimum_mag=None, maximum_mag=None, magnitude_binsize = 0.1, doplot=False):
    # a cummulative distrbution
    if minimum_mag is None:
        minimum_mag = min(M)
    if maximum_mag is None:
        maximum_mag = max(M)
        
    mags = np.arange(minimum_mag,\
                        maximum_mag+magnitude_binsize, \
                        magnitude_binsize).tolist();
    N = []
    for mbin in mags:
        n = sum(m>=mbin for m in M)
        N.append(n)
    if doplot:
        plt.semilogy(mags, N, 'o')
        plt.xlabel('Magnitude Mw')
        plt.ylabel('N')
    return (N, mags)

def calc_bvalue(mags, magbin = 0.1, Mc=None):
    if Mc is None:
        Mc = min(mags)
        M = mags
    else:
        M = [m for m in mags if m >=Mc] 
    # Aki 1965; Bender 1983; Utsu 1999
    bvalue = np.log10(np.exp(1))/(np.mean(M) - (Mc-magbin/2.0))
    return bvalue

#
def calc_avalue(bvalue, M, N):
    # N is number of events with magnitude>=M
    avalue = np.log10(N)-bvalue*M
    return avalue
   
def sample_GRdistr(bvalue, avalue, Mmin, Mmax, \
                   nevents = 1000, mbin=0.1, \
                   unround = False,
                   rseed=None):
    
    mags = np.arange(Mmin, Mmax+mbin, mbin)
    
    N = [10**(avalue-bvalue*mw) for mw in mags]
    # weights- is this necessary
    probs = [n/sum(N) for n in N]
    seed(rseed)
    sample_mag = choices(mags, probs, k=nevents)
    kmags = []
    # unround 
    if unround:
        magerr = np.random.normal(scale=0.05, size=len(sample_mag)).tolist()
        for m in sample_mag:
            kmags.append(m + np.random.normal(scale=0.05, size=1)[0])
    else:
        kmags = sample_mag    
    return kmags


def test_sample_GRdistr(bvalue, avalue, Mmin, Mmax, nevents=10000, rseed=None, doplot=True):
    #
    samp_mag = sample_GRdistr(bvalue, avalue, Mmin, Mmax,\
                     nevents = nevents, mbin=0.1, rseed=rseed)
    
    N_Mmin = 10**(avalue-bvalue*Mmin)

    eN, emag = distr_empirical(samp_mag, minimum_mag=Mmin, doplot=doplot);

    eN_Mmin = eN[emag.index(Mmin)]

    Nrat = eN_Mmin/N_Mmin

    N_model = [10**(avalue-bvalue*m)*Nrat for m in emag]

    bvalue_estimate = calc_bvalue(samp_mag, magbin = 0.1)

    N_estimated = [10**(avalue-bvalue_estimate*m)*Nrat for m in emag]
    
    if doplot:
       plt.semilogy(emag, N_model, 'r-')
       plt.semilogy(emag, N_estimated, 'b--')
       plt.xlim([Mmin, Mmax])
       plt.xlabel('Magnitude Mw')
       plt.ylabel('N')
       
    return samp_mag   
       
       
