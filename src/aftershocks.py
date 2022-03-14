
import random 
import numpy as np
import src.datetimedistr as td
import src.depthdistr as dd
import src.spatprobdistr as spd
import src.magfreqdistr as mfd
from datetime import datetime, timedelta
from random import choices, seed


# The empirical Bath's law states that the average difference in magnitude 
# between a mainshock and its largest aftershock is 1.2, 
# regardless of the mainshock magnitude.


# λ(t, m) = χ [(t + c)^−p]
def get_chi(mag):
    # -- not righly done!
    a = 0.2
    return 10**(mag*a)

def get_Omori_params(mag):
    # sample p-value
    pval = np.random.randn(1)*0.07+1
    pval = pval[0]
    
    #  c 
    # - suppose to be constant but anyways
    # - one sec to few seconds in units of days
    #c = 0.000694*random.randint(0, 9)
    c = 0.09
    # chi - scaling with mag:  χ ~ 10^(M)
    chi = get_chi(mag)    
    return(pval, chi, c)

def omori_law(mag, t, params=None):
    if params is None:
        params = get_Omori_params(mag)
   
    pval, chi, c = params
    lda = chi*((t-c)**(-pval))
    return lda


def get_aftershocks(mag, lon, lat, dep, date, doplot=False):
    #
    # Omori huh?
    params = get_Omori_params(mag)
    ndays, sdist = GardnerKnopoff_window(mag)
    t = range(1, round(ndays),1)
    n = [omori_law(mag, ti, params) for ti in t]
   
    nafshocks = round(sum(n)) 
    
    if nafshocks<1:
        return None, None, None, None, None
    
    
    # sample the magnitudes ------------
    dM_Bath = 1.2
    Mmax_afshocks = mag-dM_Bath
    bvalue = np.random.randn(1)*0.07+1.2
    avalue = np.log10(nafshocks)-bvalue*3.7  # assuming Mmin 3.7 
    # logN = a-bM
    M_afshocks = mfd.sample_GRdistr(bvalue, avalue, 3.7, Mmax_afshocks, \
                                nevents = nafshocks, mbin=0.1)   
    M_afshocks = [round(m,2) for m in M_afshocks]
   
    
    # get spatial and time bounds ---------------------
    
    # depth distribution 
    #  - aftershock flag is set True
    depth_afshocks = dd.sample_depths(nafshocks, \
                                mean_depth = dep, isaftershock=True)
    # spatial-time extents
    
    # sample the dates based on Omori
    t = range(1, round(ndays),1)
    n = [omori_law(mag, ti, params) for ti in t]
    probs = [x/sum(n) for x in n]
   
    delta_times = choices(t, probs, k=nafshocks)
    dates_afshocks = []
    for dt in delta_times:
        dates_afshocks.append(date + timedelta(days=round(dt)))
    
    # create a spatail bounds of 
    # L = fshocks_dist*2, W = (fshocks_dist*0.667)*2 in kms
    dL = sdist/111
    dW = (sdist*0.6667)/111
    
    #min_lat, max_lat = lat-dW, lat+dW
    # min_lon, max_lon = lon-dL, lon+dL
    #geobounds = (min_lon, max_lon, min_lat, max_lat)
    #lon_afshocks, lat_afshocks = spd.get_epicenters(geobounds, patch='high',\
    #                                nevents=nafshocks, \
    #                                corr_dist = round(dL),doplot=doplot)
    # get_aftershocks_epicenter(nevents, lon, lat, lx, wy, cx, cy, doplot=False):
    lon_afshocks, lat_afshocks=spd.get_aftershocks_epicenter(nafshocks, lon, lat, dL*10, dW*10, dL, dW)
    
    return M_afshocks, lon_afshocks, lat_afshocks, depth_afshocks, dates_afshocks


def get_foreshocks(mag, lon, lat, dep, date, doplot=False):
    # mag, lon, lat, dep - magnitude and hypocenter
    dM_Bath = 1.2
    M_range = 2.0 # say two units 
    # there is no general theoretical or empirical rule for 
    # total number of events either for foreschoks 
    candies = [int(x) for x in range(round(mag/2))]
    Mmax_fshocks = mag-dM_Bath

    nfshocks = random.choice(candies)
    
    if not nfshocks:
        return None, None, None, None, None
    
    # sample the magnitudes from an uniform distrbution ------------
    M_fshocks = np.random.uniform(size = nfshocks, \
                    low = Mmax_fshocks-M_range, high = Mmax_fshocks)
    
    M_fshocks = [round(m,2) for m in M_fshocks]
    
    nfshocks = len(M_fshocks)
    
    # get spatial and time bounds ---------------------
    
    # depth distribution 
    #  - aftershock flag is considered for foreshocks as well
    depth_fshocks = dd.sample_depths(nfshocks, \
                                mean_depth = dep, isaftershock=True)
    # spatial-time extents
    # we use the aftershock window, but consider much smaller domain
    ndays, sdist = GardnerKnopoff_window(mag)
    # dates
    fshocks_ndays = ndays/30
    start_date = date - timedelta(days=fshocks_ndays)
    end_date = date
   
    dates_fshocks = td.get_random_dates(nfshocks,\
                        start_date, end_date)
    #
    fshocks_dist = sdist/3
    # create a spatail bounds of 
    # L = fshocks_dist*2, W = (fshocks_dist*0.667)*2 in kms
    dL = fshocks_dist/111
    dW = (fshocks_dist*0.6667)/111
    
    min_lat, max_lat = lat-dW, lat+dW
    min_lon, max_lon = lon-dL, lon+dL
    geobounds = (min_lon, max_lon, min_lat, max_lat)
    lon_fshocks, lat_fshocks = spd.get_epicenters(geobounds, patch='high',\
                                    nevents=nfshocks, \
                                    corr_dist = round(dL),doplot=doplot)
    
    return M_fshocks, lon_fshocks, lat_fshocks, depth_fshocks, dates_fshocks


def GardnerKnopoff_window(mag):
    # Mag-Time-Distance window of Gardner and Knopoff, 1974
    # for aftershocks
    # T time is in days, L is in Km
    
    if mag>=6.5:
        log10T = 0.032*mag+2.7389
    else:
        log10T = 0.5409*mag-0.547
    
    log10L = 0.1238*mag+0.983
    
    return 10**log10T, 10**log10L


