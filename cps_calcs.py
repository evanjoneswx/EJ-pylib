import xarray as xr
import numpy as np

def b_calcs(height,center_lon,center_lat,ang):
    '''Calculates the thickness asymmetry across a TC for use within the cyclone phase space (CPS) (methodology from Hart 2003)
        Inputs: height: 3D array of geopotential height values in lat/lon
                center_lon, center_lat: Longitude and latitude of TC center
                ang: storm motion bearing
        Outputs: Thickness asymmetry value (B) for a TC for analysis within the CPS (one value in m)'''
    # constants
    d2r = np.pi / 180 # convert from degrees to radians
    h = 1.0 # h needs to be -1 if in SH

    # select the data to look at
    height_600 = height.sel(level=600)
    height_900 = height.sel(level=900)
    thickness = height_600-height_900

    # make meshgrids of the selected thickness, longitude and latitude
    thickness_subset, lon2d, lat2d = sel_subset(thickness,center_lon,center_lat)
    
    # find the bearing between the center and each lat/lon in the array
    angles_all = bearing(center_lat,center_lon,lat2d,lon2d)

    # find the 500-km radius for Z calculations and only use points inside 500-km radius
    d = haversine(center_lat, center_lon, lat2d, lon2d)
    Zl = xr.where(d < tc_rad,thickness_subset,np.nan) # set points outside 500-km radius to nans
    Zr = xr.where(d < tc_rad,thickness_subset,np.nan) # set points outside 500-km radius to nans
    # set values along great circle line to nans (since they wouldn't be either in the left or right hemisphere of the circle technically)
    Zl = xr.where(angles_all == ang,np.nan,Zl)
    Zr = xr.where(angles_all == ang,np.nan,Zr)
    # for storm motion angles in quadrants 1 & 2 (NE and SE)
    if ang >= 0 and ang < 180:
        Zl = xr.where((angles_all < ang) | (angles_all > ang+180), Zl, np.nan)
        Zr = xr.where((angles_all > ang) & (angles_all < ang + 180), Zr, np.nan)
    # for storm motion angles in quadrants 3 & 4 (NW and SW)
    elif ang >= 180 and ang < 360:
        Zl = xr.where((angles_all > ang - 180) & (angles_all < ang), Zl, np.nan)
        Zr = xr.where((angles_all > ang) | (angles_all < ang - 180), Zr, np.nan)

    # make array for calculating the weighted average
    weighted_lat = np.cos(np.deg2rad(lat2d))
    # only have weights where there are valid values
    Zr_weights = xr.where(np.isnan(Zr),np.nan,weighted_lat)
    Zl_weights = xr.where(np.isnan(Zl),np.nan,weighted_lat)
    # get the weighted value of Zr and Zl
    weighted_Zr = Zr*Zr_weights
    weighted_Zl = Zl*Zl_weights
    # calculate the weighted sum of weighted Zr and weighted Zl
    Zr_sum = weighted_Zr.sum(dim=('latitude','longitude'))
    Zl_sum = weighted_Zl.sum(dim=('latitude','longitude'))
    # sum the weights for dividing by in weighted mean for each side
    Zr_weights_sum = Zr_weights.sum(dim=('latitude','longitude'))
    Zl_weights_sum = Zl_weights.sum(dim=('latitude','longitude'))
    # calculate the weighted average of each side
    Br_weighted = Zr_sum / Zr_weights_sum
    Bl_weighted = Zl_sum / Zl_weights_sum
    # calculate B
    B = Br_weighted - Bl_weighted 

    return B


def calc_VltVut(height_all_levs,center_lon,center_lat,lev):
    '''Calculates the lower tropospheric and upper tropospheric thermal wind for a TC for use within the cyclone phase space (CPS) (methodology from Hart 2003)
        Inputs: height_all_levs: 3D array of geopotential height values in lat/lon
                center_lon, center_lat: Longitude and latitude of TC center
                lev: 1D array of level values based on height_all_levs
        Outputs: Upper tropospheric thermal wind (Vut) and lower tropospheric thermal wind (Vlt) for a TC for analysis within the CPS (one value for each in m)'''
    # constants 
    d2r = np.pi / 180 # convert from degrees to radians
    search_box = 15
    tc_rad = 500
    res = 0.25 # resolution of data
    h = 1.0 # h needs to be -1 if in SH

    # make meshgrids of the selected heights, longitude and latitude
    height_subset, lon2d, lat2d = sel_subset(height_all_levs,center_lon,center_lat)

    lon3d = np.repeat(lon2d[np.newaxis,:,:],lev.size,axis=0)
    lat3d = np.repeat(lat2d[np.newaxis,:,:],lev.size,axis=0)
    # find 500-km radius for thickness calculations
    d = haversine(center_lat, center_lon, lat3d, lon3d)
    thickness_500km = xr.where(d < tc_rad,height_subset,np.nan) # set points outside 500-km radius to nans
    dZ = thickness_500km.max(dim=['latitude','longitude']) - thickness_500km.min(dim=['latitude','longitude'])
    # select slices between 300-600 and between 600-900
    dZu = dZ.sel(level=slice(300,600))
    dZl = dZ.sel(level=(650,700,750,800,850,900))
    # take the natural log of each pressure level
    lnpu = np.log(lev.sel(level=slice(300,600)))
    lnpl = np.log(lev.sel(level=(650,700,750,800,850,900))) 
    # Compute thermal wind using linear regressions
    Vut, bupper = linear_regression(lnpu,dZu)
    Vlt, blower = linear_regression(lnpl,dZl)

    return Vut, Vlt

def weighted_mean(array):
    '''Calculates the 24-hr running mean of array of values (used here for smoothing the cyclone phase space (CPS) parameters (methodology from Hart 2003)
        Inputs: array of values to smooth over (needs to be in 6-hourly intervals for this specific smoothing function)
        Outputs: 24-hr running mean (smoothed) values of inputted array'''
    array = np.array(array)
    array_running_avg = array[0:2]
    for i in range(2,len(array)-2):
        array_5_times = array[i-2:i+3]
        array_avg = np.sum(array_5_times) / 5
        array_running_avg = np.append(array_running_avg,array_avg)

    array_end = array[len(array)-2:len(array)]
    array_running_avg = np.append(array_running_avg,array_end)
    return array_running_avg