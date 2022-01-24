import xarray as xr
import metpy.calc 
from metpy.units import units
import numpy as np

def sel_subset(array,center_lon,center_lat,search_box):
    '''Selects a subset of data relative to a central reference point and specified search area
        Inputs: array: gridded dataset to make the selection from
                center_lon, center_lat: reference center to select from
                search_box: area to select over (rectangle whose width is dimensions of search_box x search_box)
        Outputs: Array subset chosen (will be square in shape'''
    # select a subset of lats/lons within search_box / 2 of center
    lon_min = center_lon - (search_box/2)
    lon_max = center_lon + (search_box/2)
    lat_min = center_lat - (search_box/2)
    lat_max = center_lat + (search_box/2)
    if center_lon >= 360 - (search_box/2): # if there is a center_lon > 350
        remainder = 360 - center_lon 
        array_subset_1 = array.sel(longitude=slice(center_lon-(search_box/2),360),latitude=slice(center_lat+(search_box/2),center_lat-(search_box/2))) 
        array_subset_2 = array.sel(longitude=slice(0,(search_box/2) - remainder),latitude=slice(center_lat+(search_box/2),center_lat-(search_box/2))) 
        array_subset = xr.concat([array_subset_1,array_subset_2],dim='longitude')
    elif center_lon < (search_box/2): # less than 7.5deg
        remainder = search_box/2 - center_lon 
        array_subset_1 = array.sel(longitude=slice(0,center_lon+(search_box/2)),latitude=slice(center_lat+(search_box/2),center_lat-(search_box/2)))
        array_subset_2 = array.sel(longitude=slice(360-remainder,360),latitude=slice(center_lat+(search_box/2),center_lat-(search_box/2)))
        array_subset = xr.concat([array_subset_1,array_subset_2],dim='longitude')
    else:
        array_subset = array.sel(longitude=slice(lon_min,lon_max),latitude=slice(lat_max,lat_min))

    return array_subset

def temp_gradient_fctn(temp_array):
    '''Calculates the temperature gradient based on an array of gridded temperature data in lat/lon
        Inputs: temp_array: gridded dataset of temperatures
        Outputs: array of calculated temperature gradient at each grid point [K/km] '''
    lat = temp.latitude.values
    lon = temp.longitude.values
    dx, dy = metpy.calc.lat_lon_grid_deltas(lon,lat)
    dx = dx / 1000
    dy = dy / 1000

    temp_gradient = metpy.calc.gradient(temp * units.K,deltas=(dy,dx) * units.kilometer)
    temp_grad = np.sqrt((temp_gradient[0]**2) + (temp_gradient[1]**2))
    temp_grad_mag = xr.DataArray(temp_grad, dims=['latitude','longitude'], coords=[lat,lon], name='temp_grad')
    return temp_grad_mag

