import numpy as np
import xarray as xr
from pymatreader import read_mat

def mat_to_nc(mat_file_location,mat_lon,mat_lat,mat_time,mat_var_name,var_name,savename):
    '''Converts a .mat file (having dimensions of time, lat, lon here) into xarray
        Inputs: mat_file_location: path and name of mat file to open
                mat_lon, mat_lat, mat_time, mat_var_name: names of longitude, latitude, time and variable name in .mat file, respectively
                var_name: name you want to use when you create the xarray
                savename: path and filename of file to be saved
        Outputs: xarray of order time, lat, lon of original .mat file that can be saved as .nc file'''
    # open mat file for reading
    ds_mat = read_mat(mat_file_location) 

    # pull in variables as numpy arrays
    longitude = ds_mat[mat_lon]
    latitude = ds_mat[mat_lat]
    time = ds_mat[mat_time]
    var = ds_mat[mat_var_name] 

    # Format the time coordinate so that is it in the normal xarray configuration
    # note: this may have to be modified based on the time vectoring used
    # in the mat file 
    time_updated = (time * np.timedelta64(1,'h')) + np.datetime64('1900-01-01T00:00:00Z') 

    # create an xarray DataArray of the mat file 
    var_xarray = xr.DataArray(var,dims=['longitude','latitude','time'],coords=[longitude,latitude,time_updated],name=var_name)
    # transpose to the conventional time, lat, lon order
    var_xarray = F_strength.transpose('time','latitude','longitude')

    # save to netcdf, if desired
    var_xarray.to_netcdf(savename)

    return var_xarray