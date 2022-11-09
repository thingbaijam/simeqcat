import matplotlib.pyplot as plt
import numpy as np
from random import choices, seed

# 

def distr_empirical(M, Mmin=None, Mmax=None, magbin = 0.1, doplot=False):
    # a cummulative distrbution
    if Mmin is None:
        Mmin = min(M)
    if Mmax is None:
        Mmax = max(M)
        
    mags = np.arange(Mmin,Mmax+magbin, magbin).tolist()
    
    N, Mw = [],[]
    for mbin in mags:
        n = sum(m>=mbin for m in M)
       
        if n>0:
           N.append(n)
           Mw.append(mbin)
        
    if doplot:
        plt.semilogy(Mw, N, 'o')
        plt.xlabel('Magnitude Mw')
        plt.ylabel('N')
        
    return (N, Mw)


def calc_bvalue(mags, magbin = 0.1, Mc=None):
    if Mc is None:
        Mc = min(mags)
        M = mags
    else:
        M = [m for m in mags if m >=Mc] 
    # Aki 1965; Bender 1983; Utsu 1999
    bvalue = np.log10(np.exp(1))/(np.mean(M) - (Mc-(magbin/2.0)))
    return bvalue

#
def calc_avalue(bvalue, M, N):
    # N is number of events with magnitude>=M
    a = [np.log10(n)+bvalue*m for m,n in zip(M,N)]
    return np.mean(a)
 
def sample_GRdistr(bvalue, avalue, Mmin, Mmax, \
                   nevents = 1000, mbin=0.1, \
                   rseed=None):
    #
    mags = np.arange(Mmin, Mmax+mbin, mbin)
    N = [10**(avalue-bvalue*mw) for mw in mags]
    # weights- is this necessary
    probs = [n/sum(N) for n in N]
    seed(rseed)
    sample_mag = choices(mags, probs, k=nevents)
    return sample_mag


def test_sample_GRdistr(bvalue, avalue, Mmin, Mmax, nevents=10000, rseed=None, doplot=True):
    #
    # -- not quite well done!
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
       
       
