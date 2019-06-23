import numpy as np
import pygplates

def Get_FZ_Directions(X1,Y1,X2,Y2):
    long1 = np.radians(X1)
    long2 = np.radians(X2)
    lat1 = np.radians(Y1)
    lat2 = np.radians(Y2)

    bearing = np.arctan2(np.sin(long2-long1)*np.cos(lat2), np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(long2-long1))
    bearing = np.degrees(bearing)
    bearing = (bearing + 360) % 360

    return bearing

def poles_of_rotation(to_time, from_time, delta_time,rotation_model, moving_plate, fixed_plate):

    #loop through a rotation in specific time intervals to extract poles of rotation

    #create variables
    Lats = []
    Longs = []
    Angles = []
    time_change = []

    for time in np.arange(to_time,from_time,delta_time):

        to_time = time
        from_time = time+delta_time
        stage_rotation = rotation_model.get_rotation(to_time,moving_plate,from_time,fixed_plate)

        pole_lat,pole_lon,pole_angle = stage_rotation.get_lat_lon_euler_pole_and_angle_degrees()

        #to make sure that all poles are expressed in the same hemisphere
        if pole_angle < 0:
            pole_lat = -1*pole_lat
            pole_lon = pole_lon-180
            pole_angle = -1*pole_angle


        time_change.append(from_time)
        #print 'Time interval = ',time,'-',time+delta_time,', Stage Pole Lat,Lon,Angle = %f,%f,%f ' % (pole_lat,pole_lon,pole_angle)
        Lats.append(pole_lat)
        Longs.append(pole_lon)
        Angles.append(np.radians(pole_angle))

    # These next lines are necessary becuase the answers come out in the northern hemisphere,
    # need to check convention
    Longs = np.add(Longs,180.)
    Lats = np.multiply(Lats,-1)

    return Longs, Lats, Angles, time_change
