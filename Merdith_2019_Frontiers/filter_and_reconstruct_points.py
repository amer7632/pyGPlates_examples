import numpy as np
import pygplates
import itertools
import os
global_collision_parameters = (7.0, 10.0)
threshold_velocity_delta, threshold_distance_to_boundary_per_my = global_collision_parameters

def get_plate_ID_of_points(plate_partitioner, latitude, longitude):
    latLonPoints = []
    for i,j in zip(latitude,longitude):
        inputlatLon = pygplates.LatLonPoint(i, j)
        latLonPoints.append(pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon))
    plateIDs = []
    for point in latLonPoints:
        tmp_partition = plate_partitioner.partition_point(point)
        if tmp_partition == None:
            plateIDs.append(tmp_partition)
        else:
            plateIDs.append(plate_partitioner.partition_point(point).get_feature().get_reconstruction_plate_id())

    return plateIDs

def filter_points_Nones(latitude, longitude, plateIDs, zvals):

    '''
    testing new filtering points
    removes any points with None
    '''
    for ind, i in enumerate(plateIDs):
        if i is None:
            latitude[ind] = None
            longitude[ind] = None
            zvals[ind] = None

    return latitude, longitude, zvals

def recalculate_lat_lons_for_time(latitude,
                                  longitude,
                                  plateIDs,
                                  reconstruct_to_time,
                                  reconstruct_from_time,
                                  rotation_model,
                                  age_grid_mask):
    '''
    This function recalculates the latitude and longitude
    of a point given a rotation file and a to/from time
    '''

    recon_lats = []
    recon_lons = []

    shape = np.shape(latitude) #get shape of desired array

    for i,j,k in itertools.izip(latitude.data,longitude.data, plateIDs):

        for l,m,n in itertools.izip(i,j,k):
            #print l,m,n
            inputlatLon = pygplates.LatLonPoint(l,m)
            latLonPoint = pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon)

            #to time, plate id, from time
            #make_raster_at, plate id, CCD_intersection_time
            point_rotation = rotation_model.get_rotation(int(reconstruct_to_time), n, reconstruct_from_time)

            reconstructed_point = point_rotation * latLonPoint

            recon_lats.append(reconstructed_point.to_lat_lon()[0])
            recon_lons.append(reconstructed_point.to_lat_lon()[1])
    #print len(recon_lats)
    #print len(plateIDs), 'plateIDs'
    #print len(recon_lats)
    #create arrays of recon lats/lons in the same shape as our input data with the appropriate age grid mask
    masked_recon_lats = np.ma.masked_array(np.asarray(recon_lats).reshape(shape), age_grid_mask.mask)
    masked_recon_lons = np.ma.masked_array(np.asarray(recon_lons).reshape(shape), age_grid_mask.mask)

    return masked_recon_lats, masked_recon_lons

def write_xyz_file(output_filename, output_data):
    #print output_filename
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')
    output_file.close()

def create_gpml_velocity_feature(longitude_array,latitude_array):
# function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
# Long and Lat are assumed to be 1d arrays.

    multi_point = pygplates.MultiPointOnSphere(zip(latitude_array,longitude_array))

    # Create a feature containing the multipoint feature.
    # optionally, define as 'MeshNode' type, so that GPlates will recognise it as a velocity layer
    meshnode_feature = pygplates.Feature()
    meshnode_feature.set_name('Multipoint Feature')

    meshnode_feature.set_geometry(multi_point)

    output_feature_collection = pygplates.FeatureCollection(meshnode_feature)

    return output_feature_collection

def seds_filter_points_in_polygon(latitude, longitude, zvalue, reconstructed_coastlines):
    '''
    This function returns the plateID of latitude/longitude points
    at a given time from a rotation file

    It uses the topological plate boundary network to assign a plate ID.
    So it needs pyGPlates.
    '''

    #because of how the data is formatted we need two for loops to check each point
    #because we are loading in a 2d array of longitude and latitude *each*

    filtered_lats = []
    filtered_lons = []
    filtered_zvals = []

    for x,y,z in itertools.izip(latitude, longitude, zvalue):

        inputlatLon = pygplates.LatLonPoint(x,y)

        latLonPoint = pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon)
        #print latLonPoint

        isPoly = []
        #get geometry of static polygons at reconstructed time
        for coastline in reconstructed_coastlines:
            #print coastline
            coastline = pygplates.PolygonOnSphere(coastline.get_reconstructed_geometry())
            #if int(coastline.is_point_in_polygon(latLonPoint)) == 1:
                #break#print tmp
            tmp = int(coastline.is_point_in_polygon(latLonPoint))
            isPoly.append(tmp)
        if np.sum(isPoly) == 0:

            filtered_lats.append(x)
            filtered_lons.append(y)
            filtered_zvals.append(z)

    return filtered_lats, filtered_lons, filtered_zvals

def find_polygon_of_points(latitude, longitude, polygons):
    '''
    Point in polygon test
    It uses the topological plate boundary network to assign a polygon.
    So it needs pyGPlates.
    '''

    for ind, i in enumerate(latitude):

        inputlatLon = pygplates.LatLonPoint(latitude[ind],longitude[ind])
        latLonPoint = pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon)
        #print latLonPoint

        isPoly = []
        #get geometry of static polygons at reconstructed time
        for coastline in reconstructed_coastlines:
            #print coastline
            coastline = coastline.get_reconstructed_geometry()
            #if int(coastline.is_point_in_polygon(latLonPoint)) == 1:
                #break#print tmp
            tmp = int(coastline.is_point_in_polygon(latLonPoint))
            isPoly.append(tmp)
            if tmp == 1:
                latitude[ind] = None
                longitude[ind] = None
                zvals[ind] = None

    latitude = filter(None, latitude)
    longitude = filter(None, longitude)
    zvals = filter(None, zvals)

    return latitude, longitude, zvals

def filter_points_in_polygon(latitude, longitude, zvals, reconstructed_coastlines):
    '''
    This function filters out points that are overlapping continental crust
    (once overlapping we assume they are subducted).

    It uses the topological plate boundary network to assign a plate ID.
    So it needs pyGPlates.
    '''

    for ind, i in enumerate(latitude):
        if not np.isnan(i):
            inputlatLon = pygplates.LatLonPoint(latitude[ind],longitude[ind])

            latLonPoint = pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon)
            #print latLonPoint

            isPoly = []
            #get geometry of static polygons at reconstructed time
            for coastline in reconstructed_coastlines:
                #print coastline
                coastline = coastline.get_reconstructed_geometry()
                #if int(coastline.is_point_in_polygon(latLonPoint)) == 1:
                    #break#print tmp
                tmp = int(coastline.is_point_in_polygon(latLonPoint))
                isPoly.append(tmp)
                if tmp == 1:
                    latitude[ind] = None
                    longitude[ind] = None
                    zvals[ind] = None

        latitude = filter(None, latitude)
        longitude = filter(None, longitude)
        zvals = filter(None, zvals)

    return latitude, longitude, zvals

def find_polygon_of_points(latitude, longitude, polygons):
    '''
    Point in polygon test
    It uses the topological plate boundary network to assign a polygon.
    So it needs pyGPlates.
    taken from simon williams github
    '''
    polygons_of_points = []
    for ind, i in enumerate(latitude):
        #print ind,i
        inputlatLon = pygplates.LatLonPoint(latitude[ind],longitude[ind])
        latLonPoint = pygplates.convert_lat_lon_point_to_point_on_sphere(inputlatLon)

        for polygon in polygons:
            polygon = polygon.get_resolved_geometry()
            tmp = int(polygon.is_point_in_polygon(latLonPoint))
            if tmp == 1:
                polygons_of_points.append(polygon)

    return polygons_of_points

def detect_collision(timestep,
                     curr_point,
                     prev_point,
                     curr_stage_rotation,
                     prev_stage_rotation,
                     curr_topology_plate_id,
                     prev_resolved_boundary,
                     prev_topology_plate_id):
    '''
    this is taken from Simon Williams github and bastardised (slightly)
    '''

            #missing stage rotations?

    #we want to return true and false to filter points
    #so if point (a) at time 0 is at x, then at time 1, is at x + distance, that distance should
    #a reflection of normal evolution of ocean crust. if it's super big, then it's probably been subducted
    #threshold =
    #time = 100 (as in, rotating points from 101 to 100)
    prev_location_velocity = pygplates.calculate_velocities(curr_point,
                                                            prev_stage_rotation,
                                                            timestep,
                                                            pygplates.VelocityUnits.kms_per_my)[0]
    curr_location_velocity = pygplates.calculate_velocities(curr_point,
                                                            curr_stage_rotation,
                                                            timestep,
                                                            pygplates.VelocityUnits.kms_per_my)[0]

    delta_velocity = curr_location_velocity - prev_location_velocity
    delta_velocity_magnitude = delta_velocity.get_magnitude()
    # Since all feature types use the same collision parameters we can use the boundary polygon instead of iterating over its sub-segments.
    if detect_collision_using_collision_parameters(timestep,
                                                   delta_velocity_magnitude,
                                                   prev_point,
                                                   prev_resolved_boundary,
                                                   threshold_velocity_delta,
                                                   threshold_distance_to_boundary_per_my):
        return True

    return False

def detect_collision_using_collision_parameters(
                                    timestep,
                                    delta_velocity_magnitude,
                                    prev_point,
                                    prev_boundary_geometry,
                                    threshold_velocity_delta,
                                    threshold_distance_to_boundary_per_my):
        '''
        this is taken from Simon Williams github and bastardised (slightly)
        '''


        if delta_velocity_magnitude > threshold_velocity_delta:
            # Add the minimum distance threshold to the delta velocity threshold.
            # The delta velocity threshold only allows those points that are close enough to the boundary to reach
            # it given their current relative velocity.
            # The minimum distance threshold accounts for sudden changes in the shape of a plate boundary
            # which are no supposed to represent a new or shifted boundary but are just a result of the topology
            # builder/user digitising a new boundary line that differs noticeably from that of the previous time period.
            distance_threshold_radians = ((threshold_distance_to_boundary_per_my + delta_velocity_magnitude)
                * timestep / pygplates.Earth.equatorial_radius_in_kms)

            if pygplates.GeometryOnSphere.distance(prev_point,
                                                   prev_boundary_geometry,
                                                   distance_threshold_radians) is not None:
                # Detected a collision.
                return True

        return False


def exporting(x,y,z, time, savedir, POSR=True):
    '''
    this function just neatens the exporting process because it's messy
    '''
    lon = x
    lat = y
    zvals = z
    export_list_C = []
    export_list_serp = []
    export_list_total = []
    export_list_upper = []
    export_list_lower = []
    export_list_dykes = []
    export_list_gabbro = []
    export_list_SR = []
    export_list_age = []

    #0 C
    #1 serp
    #2 total thick
    #3 upper thick
    #4 lower thick
    #5 trans thick
    #6 dyke thick
    #7 gabbro thick
    #8 spreading rate
    #9 age

    for i,j,k in itertools.izip(lon,lat,zvals):
        export_list_C.append((i,j,np.float(k[0])))
        export_list_serp.append((i,j,np.float(k[1])))
        export_list_total.append((i,j,np.float(k[2])))
        export_list_upper.append((i,j,np.float(k[3])))
        export_list_lower.append((i,j,np.float(k[4])))
        export_list_dykes.append((i,j,np.float(k[5]+k[6])))#trans part of dykes?
        export_list_gabbro.append((i,j,np.float(k[7])))
        export_list_SR.append((i,j,np.float(k[8])))
        export_list_age.append((i,j,np.float(k[9]-time)))

    if POSR:
        print 'POSR'
        key='POSR'
    else:
        print 'PMSR'
        key='PMSR'
    write_xyz_file('%s%s_C_%s_Ma' % (savedir,key, time), export_list_C)
    write_xyz_file('%s%s_Serp_%s_Ma' % (savedir,key,time), export_list_serp)
    write_xyz_file('%s%s_total_%s_Ma' % (savedir,key,time), export_list_total)
    write_xyz_file('%s%s_upper_%s_Ma' % (savedir,key,time), export_list_upper)
    write_xyz_file('%s%s_lower_%s_Ma' % (savedir,key,time), export_list_lower)
    write_xyz_file('%s%s_dykes_%s_Ma' % (savedir,key,time), export_list_dykes)
    write_xyz_file('%s%s_gabbro_%s_Ma' % (savedir,key,time), export_list_gabbro)
    write_xyz_file('%s%s_SR_%s_Ma' % (savedir,key,time), export_list_SR)
    write_xyz_file('%s%s_age_%s_Ma' % (savedir,key,time), export_list_age)

    return

def write_xyz_file(output_filename, output_data):
    #print output_filename
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')
    output_file.close()

#set blockmedian tmp files
