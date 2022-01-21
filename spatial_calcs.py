import numpy as np

def haversine(lat1, lon1, lat2, lon2):  
    '''This function calculates the Haversine distance on Earth's surface (great circle distance) between two points 
        Inputs: lat1, lon1 are the reference coordinates to start with (0 360 for longitude)
                lat2, lon2 are the coordinates the distance is being calculated to. 
                Note that lat2/lon2 can be 2D meshgrids as well and would return an array of the distances
        Outputs: Haversine distance in km (either one value or array of values depending on inputs chosen'''
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