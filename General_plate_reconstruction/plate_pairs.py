import pygplates
import numpy as np

#
#these scripts were written by Simon williams
#(https://github.com/siwill22)
#and slightly adopted by me
#

def get_unique_plate_pairs_from_rotation_model(rotation_model,recon_time):
    # given a rotation model and a specifief reconstruction time, return
    # a list of the unique plate pairs

    tree = rotation_model.get_reconstruction_tree(recon_time)
    edges = tree.get_edges()
    # Get a list of plate pairs
    tree_list = []
    for edge in edges:
        if edge.get_parent_edge() is not None:
            tree_list.append((edge.get_fixed_plate_id(),edge.get_parent_edge().get_fixed_plate_id()))
            #if there's only one physical plate in the hierarchy (i.e, 0<--1<--plateID, with 0 being spin axis and
            #1 being a TPW correction for example) then it won't plot those plates
            #this if statement captures these plates
            if edge.get_parent_edge().get_fixed_plate_id() == 0:
                #print 'first', edge.get_fixed_plate_id(), edge.get_moving_plate_id()
                tree_list.append((edge.get_moving_plate_id(),edge.get_fixed_plate_id()))

    uniq_plate_pairs_from_rotations = list(set(tree_list))

    return uniq_plate_pairs_from_rotations


def get_unique_plate_ids_from_reconstructed_features(reconstructed_features):
    # given a set of reconstructed features, return a list of the unique plate ids

    feature_plate_ids = []
    for reconstructed_feature in reconstructed_features:
        feature_plate_ids.append(reconstructed_feature.get_feature().get_reconstruction_plate_id())
    unique_plate_ids = list(set(feature_plate_ids))

    return unique_plate_ids


# function to get centroid from every polygon in the reconstructed static polygons
def GetPolygonCentroid(static_polygons,plateid):
    centroid = []
    target_polygon_area = 0
    for polygon in static_polygons:
        if polygon.get_feature().get_reconstruction_plate_id()==plateid:
            if polygon.get_reconstructed_geometry() is not None:
                if polygon.get_reconstructed_geometry().get_area()>target_polygon_area:
                    centroid = polygon.get_reconstructed_geometry().get_boundary_centroid().to_lat_lon()
                    target_polygon_area = polygon.get_reconstructed_geometry().get_area()

    return centroid


# Alternatively, get centroids of topological polygons
def GetPlateCentroid(resolved_polygons,plateid):
    centroid = []
    for polygon in resolved_polygons:
        if polygon.get_feature().get_reconstruction_plate_id()==plateid:
            centroid = polygon.get_resolved_boundary().get_boundary_centroid().to_lat_lon()

    return centroid


def patch_links_between_polygon(moving_plate,uniq_rotation_pairs,uniq_plates_from_static_polygons):
    # for a given plate id, find the next highest plate id in the hierarchy
    # for which a geometry exists in the specified set of polygons
    # the first entry in the returned list is always the input moving plate
    plate_chain_list = [moving_plate]

    fixed_plate = None
    found_the_end = False
    while not found_the_end:

        # first, find the plate pair (from the rotation tree) for the input moving plate.
        # Append the fixed plate to our list
        for plate_pair in uniq_rotation_pairs:
            if plate_pair[0]==moving_plate:
                fixed_plate = plate_pair[1]
                plate_chain_list.append(fixed_plate)
                continue

        # if fixed plate is still None, we didn't find a valid plate pair - exit
        if fixed_plate is None:
            found_the_end = True
        # if fixed plate id has a valid geoemtry, we're done - exit
        elif fixed_plate in uniq_plates_from_static_polygons:
            found_the_end = True
        # if fixed plate id doesn't have a valid id, set the fixed plate to be moving plate
        # and start the loop again
        else:
            moving_plate = fixed_plate

        # Note that this is important to stop infinite loop when the fixed plate becomes zero
        # Not sure why the 'None' doesn't kick in?
        if fixed_plate==0:
            break

    return plate_chain_list


def get_plate_chains(uniq_plates_from_static_polygons,uniq_plate_pairs_from_rotations):
    # given a list of plate ids found in static polygons, and a list of plate pairs
    # from a rotation tree, finds the linkages between plate geometries, and
    # returns them as a list of lists.
    # where two plates (for which geometries exist) are linked by one or more
    # intermediate plate ids for which no geometry exists, the returned list will
    # contain three or more plate ids, where only the first and last entries
    # have valid geometries in the polygon set
    # where a plate is moving wrt to the spin axis, the list ends with zero

    chains = []
    for plate in uniq_plates_from_static_polygons:
        # call function to patch links in plate chain
        chain = patch_links_between_polygon(plate,uniq_plate_pairs_from_rotations,uniq_plates_from_static_polygons)
        # if chain is length one, means that no fixed plate could be found for
        # this moving plate (ie no point having this polygon)
        if len(chain)>1:
            #print 'no fixed plate found for moving plate id %d' % plate
            chains.append(chain)

    return chains


def create_hierarchy_features(chains,reconstructed_static_polygons,tree_features=None,valid_time=None):
    #take plate chains and static polygons, and create a set of line features that
    # join up the centroid points of polygons based on their linkage in the rotation
    # hierarchy.
    # If tree_features in given as an existing list of features, the
    # new features will be appended. Otherwise a new feature is created.
    # valid time (optional) can be given as a tuple

    if tree_features is None:
        tree_features = []

    for chain in chains:

        p0 = GetPolygonCentroid(reconstructed_static_polygons,chain[0])
        p1 = GetPolygonCentroid(reconstructed_static_polygons,chain[-1])

        if (len(p0)>0) & (len(p1)>0):  # in theory this if statement not needed now?

            feature = pygplates.Feature()
            simple_line = pygplates.PolylineOnSphere([p0,p1])
            feature.set_geometry(simple_line.to_tessellated(np.radians(1)))
            feature.set_name(str(chain))
            if valid_time is not None:
                feature.set_valid_time(valid_time[0],valid_time[1])

            tree_features.append(feature)

    return tree_features
