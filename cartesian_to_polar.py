import xarray as xr
import numpy as np

def polar_coords(var, centerlat, centerlon,current_bearing=None,name,deg_res,radius,radius_res):
    '''Converts Cartesian coordinates of latitude/longitude to polar coordinates (radius, theta (bearing))
        Note: needs Haversine and bearing formula in other file
        Inputs: var: array to be converted
                centerlat, centerlon: reference about which to convert the coordinates
                current_bearing: if rotating based off a vector, this is the angle to rotate with respect to
                name: name to save the new converted array 
                deg_res: resolution for binning in degrees
                radius: radius from the centerlat/centerlon to go out to for the polar coordinates
                radius_res: resolution for binning in km
        Outputs: xarray of variable in polar coordinates (radius, angle)'''
    var = var.to_dataset() #have to convert from DataArray to Dataset for stuff below
    # make meshgrids of lon and lat for calculations
    lon2D, lat2D = np.meshgrid(var.longitude.values, var.latitude.values)
    #now we add distance as a variable to dataset (dist only has dims [lat,lon], not lev or time)
    var['dist'] = xr.DataArray(haversine(centerlat, centerlon, lat2D, lon2D),dims=['latitude','longitude'],coords=[var.latitude,var.longitude])
    var['bearing'] = xr.DataArray(bearing(centerlat, centerlon,lat2D,lon2D),dims=['latitude','longitude'],coords=[var.latitude,var.longitude])
    
    # include these if you want to rotate the direction of the storm
    #var['bearing'] = var['bearing'] + (360-current_bearing)
    #var['bearing'] = xarray.where(var['bearing'] > 360,var['bearing'] - 360,var['bearing'])

    # convert to a pandas dataframe to make the ensuing calculations easier
    df = var.to_dataframe().reset_index()[[name,'dist','bearing']] #only need these columns
    bearing_bin, bearing_coord = np.arange(0,360+deg_res,deg_res), np.arange(deg_res,360+deg_res,deg_res)
    dist_bin, dist_coord = np.arange(0,radius+radius_res,radius_res), np.arange(radius_res,radius+radius_res,radius_res)
    df['dist_bin'] = pd.cut(df['dist'], bins=dist_bin, labels=dist_coord, include_lowest=True)
    df['bearing_bin'] = pd.cut(df['bearing'], bins=bearing_bin, labels=bearing_coord, include_lowest=True)
    df = df[[name,'dist_bin','bearing_bin']]
    f = lambda x: np.nan if x.isnull().all() else pd.Series.mode(x)[0]
    df = df.groupby(['dist_bin','bearing_bin']).agg(f).reset_index()
    df[name] = df[name].interpolate()
    d2r = np.pi / 180
    var = xarray.DataArray(np.reshape(df[name].values, (len(dist_coord),len(bearing_coord))), dims=['dist','bearing'], coords=[dist_coord, d2r*(bearing_coord)], name=name)
    var_new = xarray.DataArray(var, dims=['dist','bearing'], coords=[dist_coord, d2r*bearing_coord], name=name)
    return var_new