#
# Python function using the Google Maps Distance Matrix API
# for distance queries.
#
# Author: David Mathias
# Date: November 3, 2018
#

import googlemaps
import time

#
# Parameters:
#   source_arg: Python tuple of source locations
#               comma-separated latitudes and longitudes
#   dest_arg:   Python tuple of destination locations
#               comma-separated latitudes and longitudes
#
# Calls gmaps distance_matrix on each element of sources X dests
#
# In the data structure returned from gmaps: 
#   rows: number of sources
#   elements: number of destinations
#
# Returns:      Python list containing two lists:
#               list of distances
#               list of drive times
#
def gmap_dist(source_arg, dest_arg):

    # create lists of sources and desintations to pass to gmaps call
    sources =[]
    dests = []
    
    # get sources from the parameter tuple and
    # put in the list
    i = 0
    while i < len(source_arg):
        s = str(source_arg[i])
        s = s + ','
        s = s + str(source_arg[i+1])
        sources.append(s)
        i += 2

    # get destinations from the parameter tuple and
    # put in the list
    i = 0
    while i < len(dest_arg):
        s = str(dest_arg[i])
        s = s + ','
        s = s + str(dest_arg[i+1])
        dests.append(s)
        i += 2
    
    # create the Google maps object
    gmaps = googlemaps.Client(key='AIzaSyB15mb3E_OKrw-eXP_E8Pv1L9xpxRTJjIs')

    print("Python: about to call maps API")
    print("  sources: {}".format(sources))
    print("  dests: {}".format(dests))
    
    # call the gmaps distance matrix function
    # the returned result is a JSON object
    distance_result = gmaps.distance_matrix(sources, dests)
    
    print("Python: after call to maps API")
    print("status: {}".format(distance_result['status']))
    attempts = 1
    while (distance_result['status'] != 'OK') and (attempts < 5):
        time.sleep(1)
        distance_result = gmaps.distance_matrix(sources, dests)
        attempts += 1
        
    # This is a JSON result returned by an earlier call
    # It is included here for testing purposes to reduce calls
    # made to the API
    # distance_result = {u'status': u'OK', u'rows': [{u'elements': [{u'duration': {u'text': u'22 mins', u'value': 1304}, u'distance': {u'text': u'17.6 km', u'value': 17606}, u'status': u'OK'}]}], u'origin_addresses': [u'5284 N Innsbruck Rd, Onalaska, WI 54650, USA'], u'destination_addresses': [u'212 15th St N, La Crosse, WI 54601, USA']}
    
    # distance_result = {u'status': u'OK', u'rows': [{u'elements': [{u'duration': {u'text': u'16 mins', u'value': 949}, u'distance': {u'text': u'15.3 km', u'value': 15316}, u'status': u'OK'}, {u'duration': {u'text': u'21 mins', u'value': 1282}, u'distance': {u'text': u'17.9 km', u'value': 17894}, u'status': u'OK'}]}, {u'elements': [{u'duration': {u'text': u'16 mins', u'value': 983}, u'distance': {u'text': u'13.8 km', u'value': 13782}, u'status': u'OK'}, {u'duration': {u'text': u'6 mins', u'value': 330}, u'distance': {u'text': u'2.5 km', u'value': 2529}, u'status': u'OK'}]}], u'origin_addresses': [u'5284 N Innsbruck Rd, Onalaska, WI 54650, USA', u'212 15th St N, La Crosse, WI 54601, USA'], u'destination_addresses': [u'1452 East Ave N, Onalaska, WI 54650, USA', u'112 29th St S, La Crosse, WI 54601, USA']}

    # determine the number of sources and destinations returned 
    # from teh gmap call
    num_sources = len(distance_result['rows'])
    num_dests = len(distance_result['rows'][0]['elements'])

    # create a list to hold the distance results
    dist = [0] * (num_sources * num_dests)
    
    # create a list to hold the time results
    time = [0] * (num_sources * num_dests)

    # get the results from the returned JSON object returned by gmaps
    # and put them in the dist and time lists
    k = 0
    for i in range(num_sources):
        for j in range(num_dests):
            dist[k] = (distance_result['rows'][i]['elements'][j]['distance']['value'])
            time[k] = distance_result['rows'][i]['elements'][j]['duration']['value']
            k += 1

    # memory cleanup
    # del sources
    # del dests
    # del distance_result
    
    # return a Python list object containing
    # the list of distance and list of times
    return [dist, time]

