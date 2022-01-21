import numpy as np

def haversine(lat1, lon1, lat2, lon2):  
    '''This function calculates the Haversine distance on Earth's surface (great circle distance) between two points 
        Inputs: lat1, lon1 are the reference coordinates to start with (0 360 for longitude)
                lat2, lon2 are the coordinates the distance is being calculated to. 
                Note that lat2/lon2 can be 2D meshgrids as well and would return an array of the distances
        Outputs: Haversine distance in km (either one value or array of values depending on inputs chosen)'''
    # distance between latitudes 
    # and longitudes 
    dLat = (lat2 - lat1) * np.pi / 180.0
    dLon = (lon2 - lon1) * np.pi / 180.0
    # convert to radians 
    lat1 = (lat1) * np.pi / 180.0
    lat2 = (lat2) * np.pi / 180.0
    # apply formulae 
    a = (pow(np.sin(dLat / 2), 2) + 
            pow(np.sin(dLon / 2), 2) * 
                np.cos(lat1) * np.cos(lat2)) 
    rad = 6371
    c = 2 * np.arcsin(np.sqrt(a)) 
    dist = rad * c 
    return dist 

def hypotenuse(res,lat1,lon1):
    '''This function calculates the hypotenuse distance across a grid box based on the resolution of your dataset (most useful with reanalyses)
        Inputs: res is the resolution of the dataset in degrees (note that for rectangular resolutions like MERRA2, the function would need to be modified)
                lat1, lon1 are the reference coordinates to start with (0 360 for longitude)
                Note that lat1/lon1 can be 2D meshgrids as well and would return an array of the hypotenuses
        Outputs: Hypotenuse distance in km (either one value or array of values depending on inputs chosen)'''
    latgrid = lat1-(res/2)
    longrid = lon1-(res/2)
    dLat = (latgrid - lat1) * np.pi / 180.0
    dLon = (longrid - lon1) * np.pi / 180.0
    lat1 = (lat1) * np.pi / 180.0
    latgrid = (latgrid) * np.pi / 180.0
    a = (pow(np.sin(dLat / 2), 2) + 
            pow(np.sin(dLon / 2), 2) * 
                np.cos(lat1) * np.cos(latgrid))
    rad = 6371
    c = 2 * np.arcsin(np.sqrt(a)) 
    hypotenuse = rad * c
    return hypotenuse

def bearing(lat1,lon1,lat2,lon2):
    '''This function calculates the bearing (from true north CW) on a gridded dataset between two points
        Inputs: lat1, lon1 are the reference coordinates to start with (0 360 for longitude)
                lat2, lon2 are the coordinates the bearing is being calculated to. 
                Note that lat2/lon2 can be 2D meshgrids as well and would return an array of the bearings
        Outputs: Bearing in degrees (either one value or array of values depending on inputs chosen)'''
    # constants
    d2r = np.pi / 180 # convert from degrees to radians
    r2d = (1/d2r) # convert radians to degrees

    # convert lats/lons from degrees to radians
    lat1r = lat1*d2r
    lon1r = lon1*d2r
    lat2r = lat2*d2r
    lon2r = lon2*d2r

    # compute angle of motion based on the two lats/lons
    ang = r2d*np.arctan2(np.sin((lon2r-lon1r))*np.cos(lat2r), np.cos(lat1r)*np.sin(lat2r) - np.sin(lat1r)*np.cos(lat2r)*np.cos(lon2r-lon1r))

    # convert from -180 180 to 0 360
    ang = ang % 360
    return ang