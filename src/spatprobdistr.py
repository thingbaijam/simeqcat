import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import random as rand 

def get_aftershocks_epicenter(nevents, lon, lat, lx, wy, cx, cy, doplot=False):
    # lx wy --> length and width for ellipical 
    # correlation distance (in degrees) 
    mlon, mlat, spat_prob = spread_Gaussian(lon, lat, lx, wy, cx, cy, doplot=False)
    
    spat_prob = np.array(spat_prob)
    # Sample the locations from spatial probabilities
    # Create a flat copy of the array
    flat_spd = spat_prob.flatten()
    indexes = [k for k in range(len(flat_spd))]

    # Now, sample the flat index  
    sample_index = rand.choices(indexes, flat_spd, k = nevents)

    # retrive the matching original index for the sample index
    x_index, y_index = np.unravel_index(sample_index, spat_prob.shape)

    lons, lats = [], []

    for xi,yi in zip(x_index, y_index):
        # should add some noise
        tlon = mlon[xi][yi]+ np.random.normal(0.0, 0.03, 1)
        tlat = mlat[xi][yi]+ np.random.normal(0.0, 0.03, 1)
        lons.append(round(tlon[0],4))
        lats.append(round(tlat[0],4))

    if doplot:
        fig = plt.figure(figsize=(10, 6), dpi=80)
        ax = fig.add_subplot(111)
        ax.contourf(mlon, mlat, spat_prob, cmap='plasma')
        ax.plot(lons, lats, 'ko')
        ax.set_aspect('equal')
        ax.set_xlabel('longitude($^\circ$E)');
        ax.set_ylabel('latitude ($^\circ$N)');
    
    return lons, lats
    
    

def get_epicenters(geobounds, patch=None, nevents=10, corr_dist = 50, sampdist = 0.1, doplot=False):
    #
    # corr_dist = 50  # it works - do not set very large distance, 100 degress
    # TO DO: need to write out the random field ..
    
    lon_min, lon_max, lat_min, lat_max = geobounds
    
    # prepare the mesh
    
    mlon, mlat, sampdist = get_mesh(lon_min, lon_max, lat_min, lat_max, sampdist)
    
    # get spatial probabilties
    nl, nw = len(mlon), len(mlon[0])
    center = (round(nw/2), round(nl/2))

    # correlation_length could be a  function of magnitudes, eh!
    # normgrf_pt2patch requires corr_length to be in terms of numbers 
    
    spat_prob = normgrf_pt2patch(center, (nw, nl), \
                          correlation_length=round(corr_dist/sampdist),\
                          patch=patch)
    
    spat_prob = np.array(spat_prob)
    # Sample the locations from spatial probabilities
    # Create a flat copy of the array
    flat_spd = spat_prob.flatten()
    indexes = [k for k in range(len(flat_spd))]

    # Now, sample the flat index  
    sample_index = rand.choices(indexes, flat_spd, k = nevents)

    # retrive the matching original index for the sample index
    x_index, y_index = np.unravel_index(sample_index, spat_prob.shape)

    lons, lats = [], []

    for xi,yi in zip(x_index, y_index):
        # should add some noise
        tlon = mlon[xi][yi]+ np.random.normal(0.0, 0.03, 1)
        tlat = mlat[xi][yi]+ np.random.normal(0.0, 0.03, 1)
        lons.append(round(tlon[0],4))
        lats.append(round(tlat[0],4))

    if doplot:
        fig = plt.figure(figsize=(10, 6), dpi=80)
        ax = fig.add_subplot(111)
        ax.contourf(mlon, mlat, spat_prob, cmap='plasma')
        ax.plot(lons, lats, 'ko')
        ax.set_aspect('equal')
        ax.set_xlabel('longitude($^\circ$E)');
        ax.set_ylabel('latitude ($^\circ$N)');
    
    return lons, lats



def spread_Gaussian(lon, lat, lx, wy, cx, cy, doplot=False):
    # let see how we can get a Gaussian distrbution 
    # lx and wy radius or length from center along x and y 
    spatsamp = 0.1
    gx = np.arange(0, lx, 0.1)
    gy = np.arange(0, wy, 0.1)
    
    if len(gx)*len(gy)<100:
        gx = np.linspace(0,lx, 10)
        gy = np.linspace(0,wy, 10)
    spatsamp_x = abs(gx[0]-gx[1])
    spatsamp_y = abs(gy[0]-gy[1])
    
    
    hnx = len(gx)
    hny = len(gy)               
    x = np.arange(-hnx, hnx)
    y = np.arange(-hny, hny)
    X, Y = np.meshgrid(x, y)
    cx = cx/spatsamp_x
    cy = cy/spatsamp_y
    filter_kernel = np.exp(-(((X)/cx)**2+((Y)/cy)**2))
   
    xlon = [lon-xi for xi in gx]
    xlat = [lat-yi for yi in gy]
    for xi in gx:
        xlon.append(lon+xi)
    for yi in gy:
        xlat.append(lat+yi)
    xlon.sort()
    xlat.sort()
    glon, glat = np.meshgrid(xlon, xlat)  
                   
    if doplot:
        plt.contourf(glon, glat, filter_kernel)
        plt.axis('equal')
    return glon, glat, filter_kernel

def generate_corrnoise(dim_num_pts, correlation_length=10):

    nw, nl = dim_num_pts;
    x = np.arange(-correlation_length, correlation_length)
    y = np.arange(-correlation_length, correlation_length)

    X, Y = np.meshgrid(x, y)

    D = np.sqrt(X*X + Y*Y)
    filter_kernel = np.exp(-D**2/(2*correlation_length))
  
    # Generate n-by-n grid of spatially correlated noise
    noise = np.random.randn(nl, nw)
    corrnoise = scipy.signal.fftconvolve(noise, filter_kernel, mode='same')

    return corrnoise


def normgrf(dim_num_pts, correlation_length=10, doplot=False):
    #
    corrnoise = generate_corrnoise(dim_num_pts, \
                        correlation_length=correlation_length)

    # convert the noise into PDF
    min_noise = min([min(r) for r in corrnoise])
    t = [x-min_noise for x in corrnoise]
    sum_t = sum([sum(r) for r in t])
    spd = [x/sum_t for x in t]
    
    if doplot:
        nw, nl = dim_num_pts
        plt.contourf(np.arange(nw), np.arange(nl), spd)
        plt.axis('equal')
        
    return spd

def normgrf_pt2patch(point, dim_num_pts, correlation_length=10, doplot=False, patch=None):
    # idea is to get distribution such that:
    # - there is high values at the point if patch is high
    # - there is low values at the point if patch is low
    
    x, y = point
    converge = False
    
    if correlation_length<=1:
        correlation_length = 1
    
    for i in range(10000):
        spd = normgrf(dim_num_pts, correlation_length)
        max_pd = max([max(r) for r in spd])
        if patch is None:
            break
        if spd[y][x]>= max_pd*0.667:
            converge = True
            break
            
    if patch =='low':
        spd = [1-x for x in pd]
    if not converge:
        print('*** not converged!!')
    if doplot:
        nw, nl = dim_num_pts
        plt.contourf(np.arange(nw), np.arange(nl), spd)
        plt.axis('equal')
        
    return spd


def get_mesh(lon_min, lon_max, lat_min, lat_max, sampdist = 0.1):
      
    # sampdist - spatial sampling distance
    #      - do not make this too small
    #      -  sampdist is override if number of grid points < 100
    
    y = np.arange(lat_min, lat_max, sampdist)
    x = np.arange(lon_min, lon_max, sampdist)
    
    ngridpts = len(y)*len(x)
    if ngridpts<100:
        y = np.linspace(lat_min, lat_max,10)
        x = np.linspace(lon_min, lon_max,10)
        sampdist = abs(x[0]-x[1])
        
    lon, lat= np.meshgrid(x, y)
    return lon, lat, sampdist

