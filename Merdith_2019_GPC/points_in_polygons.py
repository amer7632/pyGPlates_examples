# -*- coding: utf-8 -*-

"""
    Copyright (C) 2017 The University of Sydney, Australia

    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


#####################################################################
# Efficient point-in-polygon testing when there are many            #
# relatively uniformly spaced points to be tested against polygons. #
#####################################################################
#
#
# For example, to find the plate ID of the polygon containing each point in a sequence of points:
#
#
#    import points_in_polygons
#
#    # A list of 'pygplates.PointOnSphere' points.
#    points = [...]
#
#    # Some polygon features (eg, coastlines).
#    polygon_feature_collection = pygplates.FeatureCollection('polygons.gpml')
#
#    # Extract the polygons from the features.
#    polygons = []
#    polygon_features = []
#    for polygon_feature in polygon_feature_collection:
#        polygons.append(polygon_feature.get_geometry())
#        polygon_features.append(polygon_feature)
#
#    # Find the polygon (features) containing the points.
#    polygon_features_containing_points = points_in_polygons.find_polygons(points, polygons, polygon_features)
#
#    # Assign a plate ID to each point (or 0 if point outside all polygons).
#    plate_ids = []
#    for polygon_feature in polygon_features_containing_points:
#        if polygon_feature is not None:
#            plate_id = polygon_feature.get_reconstruction_plate_id()
#        else:
#            plate_id = 0
#
#        plate_ids.append(plate_id)
#
#####################################################################


from __future__ import print_function
import math
import sys
sys.path.insert(1, '/Users/andrew/Documents/python/pygplates_rev28_python37_MacOS64')
import pygplates
import sys
import numpy as np



def find_polygons(
        points,
        polygons,
        polygon_proxies = None,
        subdivision_depth = 4):
    """
    Efficient point-in-polygon testing when there are many relatively uniformly spaced points to be tested against polygons.

    points: a sequence of 'pygplates.PointOnSphere'.

    polygons: a sequence of 'pygplates.PolygonOnSphere'.

    polygon_proxies: Optional sequence of objects associated with 'polygons'.
                     If not specified then the proxies default to the polygons themselves.
                     These can be any object (such as the 'pygplates.Feature' that the polygon came from).

    subdivision_depth: The depth of the lat/lon quad tree used to speed up point-in-polygon queries.
                       The lat/lon width of a leaf quad tree node is (90 / (2^subdivision_depth)) degrees.
                       Generally the denser the 'points' the larger the depth should be.
                       Setting this value too high causes unnecessary time to be spent generating a deep quad tree.
                       Setting this value too low reduces the culling efficiency of the quad tree.
                       However a value of 4 seems to work quite well for a uniform lat/lon spacing of 'points' of 1 degree and below
                       without the cost of generating a deep quad tree.
                       So most of the time the subdivision depth can be left at its default value.

    Returns: A list of polygon proxies associated with 'points'.
             The length of the list matches the length of 'points'.
             For each point in 'points', if the point is contained by a polygon then that polygon's proxy
             is stored (otherwise None is stored) at the same index (as the point) in the returned list.

    Raises ValueError if the lengths of 'polygons' and 'polygon_proxies' (if specified) do not match.
    """

    quad_tree = QuadTree(points, subdivision_depth)
    return quad_tree.find_polygons(polygons, polygon_proxies)


##################
# Implementation #
##################


class QuadTree(object):

    def __init__(self, points, subdivision_depth):

        if subdivision_depth < 0:
            raise ValueError('Subdivision depth must be a non-negative value.')
        elif subdivision_depth > 100:
            raise ValueError('Subdivision depth is too large (should be 100 or less).')

        # Keep a copy of the points.
        self.points = points[:]

        # Each root quad tree node is quadrant of the globe (square in lat/lon space of size 90 x 90 degrees).
        # So there are 8 of them.
        # We'll only create them as needed.
        self.root_nodes = [None] * 8

        # Place each point in a quad tree leaf node.
        for point_index, point in enumerate(points):
            point_lat, point_lon = point.to_lat_lon()

            # Get root node that current point is in.
            if point_lat < 0:
                root_node_lat_index = 0
            else:
                root_node_lat_index = 1

            if point_lon < 0:
                if point_lon < -90:
                    root_node_lon_index = 0
                else:
                    root_node_lon_index = 1
            else:
                if point_lon < 90:
                    root_node_lon_index = 2
                else:
                    root_node_lon_index = 3

            root_node_index = 4 * root_node_lat_index + root_node_lon_index
            root_node = self.root_nodes[root_node_index]

            is_north_hemisphere = (root_node_lat_index == 1)
            root_node_half_width_degrees = 45.0
            root_node_centre_lon = -180 + 90 * root_node_lon_index + root_node_half_width_degrees
            root_node_centre_lat = -90 + 90 * root_node_lat_index + root_node_half_width_degrees

            # Create root node if first time visited.
            if root_node is None:
                root_node = self._create_node(
                        root_node_centre_lon,
                        root_node_centre_lat,
                        root_node_half_width_degrees,
                        is_north_hemisphere)
                self.root_nodes[root_node_index] = root_node

            # Iterate through the subdivision levels and place current point in the correct quad tree leaf node.
            node = root_node
            node_half_width_degrees = root_node_half_width_degrees
            node_centre_lon = root_node_centre_lon
            node_centre_lat = root_node_centre_lat
            for level in np.arange(0, subdivision_depth):

                child_node_half_width_degrees = node_half_width_degrees / 2.0

                if point_lat < node_centre_lat:
                    child_node_lat_offset = 0
                    child_node_centre_lat = node_centre_lat - child_node_half_width_degrees
                else:
                    child_node_lat_offset = 1
                    child_node_centre_lat = node_centre_lat + child_node_half_width_degrees

                if point_lon < node_centre_lon:
                    child_node_lon_offset = 0
                    child_node_centre_lon = node_centre_lon - child_node_half_width_degrees
                else:
                    child_node_lon_offset = 1
                    child_node_centre_lon = node_centre_lon + child_node_half_width_degrees

                # The current node is an internal node (because it will have child nodes).
                # Create a list of child nodes if first time visiting node.
                if node.child_nodes is None:
                    # Only create each child node as needed.
                    node.child_nodes = [None] * 4

                child_node_index = 2 * child_node_lat_offset + child_node_lon_offset
                child_node = node.child_nodes[child_node_index]

                if child_node is None:
                    child_node = self._create_node(
                            child_node_centre_lon,
                            child_node_centre_lat,
                            child_node_half_width_degrees,
                            is_north_hemisphere)
                    node.child_nodes[child_node_index] = child_node

                # Child node becomes parent node in next iteration.
                node = child_node
                node_half_width_degrees = child_node_half_width_degrees
                node_centre_lon = child_node_centre_lon
                node_centre_lat = child_node_centre_lat

            # Reached leaf node (end of subdivision).
            # Create a list of point indices if first time visiting node.
            if node.point_indices is None:
                node.point_indices = []

            # Add the current point (index) to the leaf node.
            node.point_indices.append(point_index)


    def find_polygons(self, polygons, polygon_proxies=None):

        # Use the polygons as proxies if no proxies have been specified.
        if polygon_proxies is None:
            polygon_proxies = polygons

        if len(polygons) != len(polygon_proxies):
            raise ValueError('Number of polygons must match number of proxies.')

        # Sort the polygons from largest to smallest area.
        # This makes searching for points/geometries more efficient.
        #
        # 'polygons_and_proxies' is a list of 2-tuples (polygon, polygon_proxy).
        polygons_and_proxies = sorted(
                ((polygons[index], polygon_proxies[index]) for index in np.arange(len(polygons))),
                key=lambda polygon_and_proxy: polygon_and_proxy[0].get_area(),
                reverse=True)

        # By default all points are outside all polygons.
        # If any are found to be inside then we'll set the relevant polygon proxy.
        polygon_proxies_containing_points = [None] * len(self.points)

        # Use a quad tree for efficiency - enables us to cull large groups of points that are either
        # outside all polygons or inside a polygon (avoids point-in-polygon tests for these points).
        for root_node in self.root_nodes:
            if root_node is not None:
                self._visit_node(root_node, polygons_and_proxies, polygon_proxies_containing_points)

        return polygon_proxies_containing_points


    def _create_node(self, node_centre_lon, node_centre_lat, node_half_width_degrees, is_north_hemisphere):

        # Create the points of the polygon bounding the current quad tree node.
        bounding_polygon_points = []

        left_lon = node_centre_lon - node_half_width_degrees
        right_lon = node_centre_lon + node_half_width_degrees
        bottom_lat = node_centre_lat - node_half_width_degrees
        top_lat = node_centre_lat + node_half_width_degrees

        # Northern and southern hemispheres handled separately.
        if is_north_hemisphere:
            # Northern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere([(0, right_lon), (90, right_lon)])

            # Midpoint of small circle arc bounding the bottom of quad tree node.
            bottom_mid_point = pygplates.PointOnSphere(bottom_lat, 0.5 * (left_lon + right_lon))

            # Find the great circle (rotation) that passes through the bottom midpoint (and is oriented towards North pole).
            bottom_great_circle_rotation_axis = pygplates.Vector3D.cross(
                    bottom_mid_point.to_xyz(),
                    pygplates.Vector3D.cross(pygplates.PointOnSphere.north_pole.to_xyz(), bottom_mid_point.to_xyz())
                            ).to_normalised()
            bottom_great_circle_rotation = pygplates.FiniteRotation(bottom_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi)

            # Intersect great circle bottom boundary with left and right boundaries to find bottom-left and bottom-right points.
            # The bottom boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our bottom boundary otherwise it will go inside the quad tree node.
            bottom_boundary = pygplates.PolylineOnSphere(
                    [bottom_great_circle_rotation * bottom_mid_point, bottom_mid_point, bottom_great_circle_rotation.get_inverse() * bottom_mid_point])
            _, _, bottom_left_point = pygplates.GeometryOnSphere.distance(bottom_boundary, left_boundary, return_closest_positions = True)
            _, _, bottom_right_point = pygplates.GeometryOnSphere.distance(bottom_boundary, right_boundary, return_closest_positions = True)

            bounding_polygon_points.append(bottom_left_point)
            bounding_polygon_points.append(bottom_right_point)

            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, right_lon))
            bounding_polygon_points.append(pygplates.PointOnSphere(top_lat, left_lon))
        else:
            # Southern hemisphere.
            left_boundary = pygplates.PolylineOnSphere([(0, left_lon), (-90, left_lon)])
            right_boundary = pygplates.PolylineOnSphere([(0, right_lon), (-90, right_lon)])

            # Midpoint of small circle arc bounding the top of quad tree node.
            top_mid_point = pygplates.PointOnSphere(top_lat, 0.5 * (left_lon + right_lon))

            # Find the great circle (rotation) that passes through the top midpoint (and is oriented towards North pole).
            top_great_circle_rotation_axis = pygplates.Vector3D.cross(
                    top_mid_point.to_xyz(),
                    pygplates.Vector3D.cross(pygplates.PointOnSphere.north_pole.to_xyz(), top_mid_point.to_xyz())
                            ).to_normalised()
            top_great_circle_rotation = pygplates.FiniteRotation(top_great_circle_rotation_axis.to_xyz(), 0.5 * math.pi)

            # Intersect great circle top boundary with left and right boundaries to find top-left and top-right points.
            # The top boundary is actually a small circle (due to lat/lon grid), but since we need to use *great* circle arcs
            # in our geometries we need to be a bit loose with our top boundary otherwise it will go inside the quad tree node.
            top_boundary = pygplates.PolylineOnSphere(
                    [top_great_circle_rotation * top_mid_point, top_mid_point, top_great_circle_rotation.get_inverse() * top_mid_point])
            _, _, top_left_point = pygplates.GeometryOnSphere.distance(top_boundary, left_boundary, return_closest_positions = True)
            _, _, top_right_point = pygplates.GeometryOnSphere.distance(top_boundary, right_boundary, return_closest_positions = True)

            bounding_polygon_points.append(top_left_point)
            bounding_polygon_points.append(top_right_point)

            bounding_polygon_points.append(pygplates.PointOnSphere(bottom_lat, right_lon))
            bounding_polygon_points.append(pygplates.PointOnSphere(bottom_lat, left_lon))

        bounding_polygon = pygplates.PolygonOnSphere(bounding_polygon_points)

        return QuadTreeNode(bounding_polygon)


    def _visit_node(self, node, parent_overlapping_polygons_and_proxies, polygon_proxies_containing_points):
        # See if the current quad tree node's bounding polygon overlaps any polygons.
        overlapping_polygons_and_proxies = []
        for polygon, polygon_proxy in parent_overlapping_polygons_and_proxies:

            # See if quad tree node and current polygon overlap.
            if pygplates.GeometryOnSphere.distance(
                    node.bounding_polygon,
                    polygon,
                    1e-4, # Arbitrarily small threshold for efficiency since only interested in zero distance (intersection).
                    geometry1_is_solid = True,
                    geometry2_is_solid = True) == 0:

                overlapping_polygons_and_proxies.append((polygon, polygon_proxy))

                # See if quad tree node is contained completely inside polygon.
                # We test this by only considering the quad tree node polygon as solid (the polygon is an outline).
                if pygplates.GeometryOnSphere.distance(
                        node.bounding_polygon,
                        polygon,
                        1e-4, # Arbitrarily small threshold for efficiency since only interested in zero distance (intersection).
                        geometry1_is_solid = True) != 0:

                    # Recursively fill the entire quad sub-tree as inside current polygon.
                    self._fill_node_inside_polygon(node, polygon_proxy, polygon_proxies_containing_points)
                    return

        # If quad tree is outside all polygons then nothing left to do since all points are marked as outside by default.
        if not overlapping_polygons_and_proxies:
            return

        # Visit child nodes (if internal node) or test each point (if leaf node).
        if node.child_nodes:
            for child_node in node.child_nodes:
                if child_node is not None:
                    self._visit_node(child_node, overlapping_polygons_and_proxies, polygon_proxies_containing_points)
        else:
            for point_index in node.point_indices:
                point = self.points[point_index]
                for polygon, polygon_proxy in overlapping_polygons_and_proxies:
                    if polygon.is_point_in_polygon(point):
                        # Point is inside a polygon.
                        polygon_proxies_containing_points[point_index] = polygon_proxy
                        break


    def _fill_node_inside_polygon(self, node, polygon_proxy, polygon_proxies_containing_points):
        if node.child_nodes:
            for child_node in node.child_nodes:
                if child_node is not None:
                    self._fill_node_inside_polygon(child_node, polygon_proxy, polygon_proxies_containing_points)
        else:
            for point_index in node.point_indices:
                # Point is inside a polygon.
                polygon_proxies_containing_points[point_index] = polygon_proxy


class QuadTreeNode(object):
    def __init__(self, bounding_polygon):
        self.bounding_polygon = bounding_polygon
        # If an internal quad tree node then 'child_nodes' will be a list of
        # 4 child quad tree nodes and 'point_indices' will be None.
        # Otherwise quad tree node is a leaf node where 'point_indices' is a list of points and
        # 'child_nodes' will be None.
        self.child_nodes = None
        self.point_indices = None


if __name__ == '__main__':

    #
    # Some testing/example code.
    #

    import time


    print('Loading coastline polygons and rotation model...')
    coastline_features = pygplates.FeatureCollection('../../../../sample_data/2.0/SampleData/FeatureCollections/Coastlines/Matthews_etal_GPC_2016_Coastlines.gpmlz')
    rotation_model = pygplates.RotationModel('../../../../sample_data/2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot')

    print('Reconstructing coastline polygons...')
    reconstruction_time = 0
    coastline_reconstructed_feature_geometries = []
    pygplates.reconstruct(coastline_features, rotation_model, coastline_reconstructed_feature_geometries, reconstruction_time)

    polygons = []
    polygon_features = []
    for reconstructed_feature_geometry in coastline_reconstructed_feature_geometries:
        polygons.append(reconstructed_feature_geometry.get_reconstructed_geometry())
        polygon_features.append(reconstructed_feature_geometry.get_feature())

    # Create uniform lat/lon distribution of points.
    print('Creating lat/lon grid of points...')
    num_latitudes = 180
    num_longitudes = 360
    lat_grid_spacing_degrees = 180.0 / num_latitudes
    lon_grid_spacing_degrees = 360.0 / num_longitudes

    points = []
    for lat_index in np.arange(num_latitudes):
        # The 0.5 puts the point in the centre of the grid pixel.
        # This also avoids sampling right on the poles.
        lat = -90 + (lat_index + 0.5) * lat_grid_spacing_degrees

        for lon_index in np.arange(num_longitudes):
            # The 0.5 puts the point in the centre of the grid pixel.
            # This also avoids sampling right on the dateline where there might be
            # age grid or static polygon artifacts.
            lon = -180 + (lon_index + 0.5) * lon_grid_spacing_degrees

            point = pygplates.PointOnSphere(lat, lon)
            points.append(point)

    print('Finding polygons containing points...')
    time_begin = time.clock()

    if True:
        #
        # The fast way (about 3 seconds).
        #
        polygon_features_containing_points = find_polygons(points, polygons, polygon_features)
    else:
        #
        # The slow way (about 200 seconds).
        #
        # Similar to 'find_polygons()' except without using a quad tree.
        #
        polygons_and_features = sorted(
                ((polygons[index], polygon_features[index]) for index in np.arange(len(polygons))),
                key=lambda polygon_and_feature: polygon_and_feature[0].get_area(),
                reverse=True)
        polygon_features_containing_points = [None] * len(points)
        for point_index, point in enumerate(points):
            for polygon, polygon_feature in polygons_and_features:
                if polygon.is_point_in_polygon(point):
                    polygon_features_containing_points[point_index] = polygon_feature
                    break

    time_end = time.clock()
    print('  {0} seconds'.format(time_end - time_begin))

    print('Associate each point with a polygon plate ID...')

    # Group points inside each polygon so can create one multi-point per polygon.
    polygon_feature_to_points_mapping = {}
    for point_index, polygon_feature in enumerate(polygon_features_containing_points):
        points_in_polygon = polygon_feature_to_points_mapping.setdefault(polygon_feature, [])
        points_in_polygon.append(points[point_index])

    # Create multi-point features.
    multi_point_features = []
    for polygon_feature, points_in_polygon in polygon_feature_to_points_mapping.iteritems():
        multi_point_feature = pygplates.Feature()
        multi_point_feature.set_geometry(
                pygplates.MultiPointOnSphere(points_in_polygon))

        # If points contained by any polygon then assign its plate ID, otherwise no plate ID assigned.
        if polygon_feature is not None:
            begin_time, end_time = polygon_feature.get_valid_time()
            multi_point_feature.set_valid_time(begin_time, end_time)

            multi_point_feature.set_reconstruction_plate_id(
                    polygon_feature.get_reconstruction_plate_id())
        else:
            multi_point_feature.set_valid_time(
                    pygplates.GeoTimeInstant.create_distant_past(),
                    pygplates.GeoTimeInstant.create_distant_future())

        multi_point_features.append(multi_point_feature)

    print('Writing points feature collection...')
    pygplates.FeatureCollection(multi_point_features).write('multi_point_features.gpml')
