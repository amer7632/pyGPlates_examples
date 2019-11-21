#this script was written by Suzanne Atkins (March, 2019), Suzanne.atkins@ens.fr (www.github.com/seatkins)
#parts of the script were adapted from work by Simon Williams, John Cannon and Nicky Wright
#(see https://github.com/siwill22/agegrid-0.1).

import pygplates as pygplates
import numpy as np
import matplotlib.pyplot as plt

def spread_rate(boundary_section,rotation_model,time):
    '''
    Each mid-ocean ridge constructed in GPlates contains segments (defined as
    a straight line between two digital points in GPLates) of spreading
    (i.e. ridges) and segments of oceanic transform faults (OTFs) that connect these
    segments. This script takes a boundary section of a plate model and loops through
    the segments to separate the spreading from fault segments based on the angle
    to the stage pole. It requires:

    boundary_section: a boundary section of a topological plate model
    rotation_model: the associated rotation file with the plate model
    time: the time of interest
    '''

    velocity_seg = [] # velocity x segment length
    velocity_ret = [] # velocity returned
    velocity_pl = [] # velocity per length
    seg_len_ret = [] # all segment lengths
    seg_plate_id = [] #segment plate IDs
    vel_mm = np.zeros([2]) #velocity in mm
    vel_mm[0] = np.inf


    shd_len=0

    #start loop through segments in plate boundary network
    for shared_sub_segment in boundary_section.get_shared_sub_segments():

        #set variables
        left_plate_id = shared_sub_segment.get_feature().get_left_plate()
        right_plate_id = shared_sub_segment.get_feature().get_right_plate()
        plate_id = shared_sub_segment.get_feature().get_reconstruction_plate_id()
        stage_rotation = rotation_model.get_rotation(time +1 , left_plate_id,time, right_plate_id, 0)
        finite_rotation = rotation_model.get_rotation(time, right_plate_id, 0, 0)
        stage_pole, stage_angle_radians = stage_rotation.get_euler_pole_and_angle()
        stage_pole_recon = finite_rotation * stage_pole
        shd_len = shared_sub_segment.get_geometry().get_arc_length()


        for segment in shared_sub_segment.get_geometry().get_segments():
            split_geometry = (segment.get_start_point(),segment.get_end_point())
            if segment.is_zero_length(): # check to see if segment is zero length
                continue

            # Get the point in the middle of the segment and its tangential direction.
            segment_midpoint = segment.get_arc_point(0.5)
            segment_direction_at_midpoint = segment.get_arc_direction(0.5)

            # Get the direction from the segment midpoint to the stage pole.
            # This is the tangential direction at the start of an arc from the segment
            # midpoint to the stage pole (the zero parameter indicates the arc start
            # point which is the segment midpoint).
            segment_to_stage_pole_direction = pygplates.GreatCircleArc(
                    segment_midpoint, stage_pole_recon).get_arc_direction(0)

            # The angle that the segment deviates from the stage pole direction.
            deviation_of_segment_direction_from_stage_pole = \
            pygplates.Vector3D.angle_between(segment_direction_at_midpoint,\
                                             segment_to_stage_pole_direction)

            # Change the angle to be between 0 and 90
            deviation_of_segment_direction_from_stage_pole_mod = \
            90-np.abs(90-np.remainder(np.degrees(deviation_of_segment_direction_from_stage_pole),\
                                      180))

            if deviation_of_segment_direction_from_stage_pole_mod <= 70:
                # make file for plotting on a map
                seg_plate_id.append([left_plate_id,right_plate_id])
                seg_len = segment.get_arc_length()
                reconstructed_points = []

                for geometry in split_geometry:
                    for point in geometry.get_points():
                        #print point
                        point_X = point.to_lat_lon()[1]
                        point_Y = point.to_lat_lon()[0]
                        velocity_point = ((point_Y,point_X))

                    reconstructed_points.append(velocity_point)

                #get_velocity of left and right divergent plates at mid ocean ridges
                #if no left or right plate id, get velocity of plate id (i.e. velocity relative to spin axis)
                if left_plate_id == 0 or right_plate_id == 0:

                    equivalent_stage_rotation = rotation_model.get_rotation(time, plate_id, time+1)

                    velocity = pygplates.calculate_velocities(reconstructed_points,equivalent_stage_rotation,1,\
                                                                pygplates.VelocityUnits.cms_per_yr)

                    velocity_magnitude_azimuth = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(\
                    reconstructed_points,velocity)

                else:
                    #if we do have relative plate motions we get the relative velocity
                    #(to_time, moving_plate, from_time, fixed_plate)
                    velocity = pygplates.calculate_velocities(reconstructed_points,stage_rotation,1,\
                                                                   pygplates.VelocityUnits.cms_per_yr)

                    velocity_magnitude_azimuth = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(\
                    reconstructed_points,velocity)

                seg_len_ret.append((seg_len*pygplates.Earth.mean_radius_in_kms))
                #magnitude only, times segment length so that mean can be found
                velocity_seg.append((velocity_magnitude_azimuth[0][0]*seg_len*pygplates.Earth.mean_radius_in_kms))
                velocity_ret.append((velocity_magnitude_azimuth[0][0]))
                velocity_pl.append((velocity_magnitude_azimuth[0][0])/(seg_len*pygplates.Earth.mean_radius_in_kms))
                if velocity_magnitude_azimuth[0][0] > vel_mm[1]:
                    vel_mm[1] = velocity_magnitude_azimuth[0][0]
                elif velocity_magnitude_azimuth[0][0] < vel_mm[0]:
                    vel_mm[0] = velocity_magnitude_azimuth[0][0]

    return velocity_ret,velocity_pl,velocity_seg,vel_mm,shd_len*pygplates.Earth.mean_radius_in_kms,seg_len_ret,seg_plate_id
