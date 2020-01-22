import numpy as np
import sys
sys.path.insert(0, '/Applications/GPlates-2.1.0/pygplates_rev18_python27_MacOS64/')
import pygplates as gplates

#import pandas as pd
from datetime import datetime
from scipy import interpolate
from collections import defaultdict
import dill
import pickle
import MOR_characterisation, filter_and_reconstruct_points
import os
import time as tme
import itertools


spread_dist = pickle.load(open('/Users/Andrew/Dropbox/1Mysteps_02042019/base_files_scripts/Distribution_for_random_spreading_Pacific100My.p.py','rb'))
loaddir = '/Users/Andrew/Dropbox/1Mysteps_02042019/base_files_scripts'
savedir = '/Users/Andrew/Dropbox/1Mysteps_02042019'


resolution = 1



#select random pacific spreading rate
def random_spreading(samples,distribution):
    y = np.random.randint(1001,size=samples)
    return distribution[y]

def crust_characterisation(start, stop, step, df):
  #%%time
  #this cell performs the calculation at each time step and then stores the results (along with other data)
  #needed for plotting and visualisation in a dictionary of dictionaries

  C_storage = defaultdict(lambda: defaultdict(list))

  # create (alot of)empty arrays to write to for total length and lengths of appropriate spreading rates

  ridge_length = []
  transform_length = []
  subduction_length = []
  transformgpml_length = []

  # Set up time array
  min_time = start
  max_time = stop
  time_step = step
  times = np.arange(min_time,max_time + time_step,time_step)

  # deviation angle (in degrees) to split transforms vs. ridges from
  deviation_angle= 70
  #time resolution to extract data
  time_resolution = 1

  SR = 'Spreading Rate'
  per = 'Peridotite'
  DS = 'Degree of Serpentinisation'
  thick = 'Thickness'
  vertical_area = 'Vertical Area'
  C_total = 'C-prod'
  C_total_vertical = 'C-prod_vertical'
  C_volcanic = 'C-Volcanic'
  C_serpentinites = 'C_Serpentinites'
  serpentinites_thick = 'Serpentinites_raw'
  serpentinites_vertical = 'Serpentinites_raw_vertical'
  length = 'Boundary-Length'
  boundary = 'Boundary-Type'
  StartLat = 'Start-Lat-Point'
  StartLon = 'Start-Lon-Point'
  EndLat = 'End-Lat-Point'
  EndLon = 'End-Lon-Point'
  RightPlateID = 'RightPlate'
  LeftPlateID = 'LeftPlate'
  Index = 'Index'

  panthalassa_plate_ids = ['902',
                                       '919',
                                       '926']
  sample=1
  #inititate loop to extract data
  for time in times:
      print time,  'Ma'
      #we only need to cut up spreading ridges by velocity
      subset1 = df[(df['Time_Ma']>=time)
                  & (df['Time_Ma']<(time+time_resolution))
                  & (df['FeatureType']=='gpml:MidOceanRidge')
                  & (np.abs(df['Deviation_mod_deg'])<=deviation_angle)]
      #print len(subset1)

      #create temporary storage for some values
      C_storage[time][SR] = []
      C_storage[time][per] = []
      C_storage[time][DS] = []
      C_storage[time][thick] = []
      C_storage[time][vertical_area] = []
      C_storage[time][C_total] = []
      C_storage[time][C_total_vertical] = []
      C_storage[time][C_volcanic] = []
      C_storage[time][C_serpentinites] = []
      C_storage[time][serpentinites_thick] = []
      C_storage[time][serpentinites_vertical] = []
      C_storage[time][length] = []
      C_storage[time][boundary] = []
      C_storage[time][StartLat] = []
      C_storage[time][StartLon] = []
      C_storage[time][EndLat] = []
      C_storage[time][EndLon] = []
      C_storage[time][RightPlateID] = []
      C_storage[time][LeftPlateID] = []
      C_storage[time][Index] = []

      for index, row in subset1.iterrows():


          #returns full spreading rate/velocity in cm/year, times by 10 to convert to km/Ma
          velocity = row.Plate_Velocity*10

          #block out if not using plate model
          if time > 180:
              new_conjugate_plates = (row.RightPlate,row.LeftPlate)
              if all(x in panthalassa_plate_ids for x in new_conjugate_plates) == True:
                  velocity = random_spreading(sample,spread_dist)*10
                  #print 'here', velocity

          #calculate the mean bottom water temperature for the first 20 Ma existence of a parcel of ocean crust
          temp_tmp = []
          for i in np.arange(time,time-21,-1):
              if i < 0:
                  break
              temp_tmp.append(bottom_water_temperature_curve[i])
          bottom_water_temperature = np.mean(temp_tmp)

          per = MOR_characterisation.SR_and_peridotite(velocity)
          volcanic_percent = 100 - per
          asymmetry_factor = 1
          calc_length = row.Length_km

          #calculate variables for mid ocean ridge segments

          calc_length = row.Length_km
          width = velocity
          if velocity <= 40:
              asymmetry_factor = 1
              thickness = MOR_characterisation.SR_and_thickness_slow_ultraslow(velocity)
              tmp_DS = MOR_characterisation.SR_and_dsSurf_slow_ultraslow(velocity, thickness) #DS and thickness is called within this function
              #print DS
              carbon_max, CO2_gabbro = MOR_characterisation.carbon_content_slow_ultraslow(velocity, per)
              upper_volc_CO2 = lower_volc_CO2 = transition_CO2 = sheeted_dykes_CO2 = 0
              upper_volc_thickness = lower_volc_thickness = transition_thickness = sheeted_dykes_thickness = 0
              gabbros_thickness = thickness * volcanic_percent/100
              if CO2_gabbro < 0:
                  CO2_gabbro = 0
              CO2_volcanic = gabbros_thickness * CO2_gabbro

          else:
              carbon_max = 0
              tmp_DS = MOR_characterisation.SR_and_dsSurf_inter_fast(velocity) #DS and thickness is called within this function

              thickness, upper_volc_thickness, lower_volc_thickness,transition_thickness,sheeted_dykes_thickness,gabbros_thickness = MOR_characterisation.SR_and_thickness_inter_fast(
                                                                                                                      tmp_DS,
                                                                                                                      volcanic_percent)
              upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2, CO2_gabbro = MOR_characterisation.carbon_content_inter_fast(velocity,
                                                                                                             per,
                                                                                                             bottom_water_temperature)
              if CO2_gabbro < 0:
                  CO2_gabbro = 0
              #gives us total CO2 storage in volcanics as a fraction (NOT A PERCENT)
              CO2_volcanic =  (upper_volc_thickness * upper_volc_CO2 + \
                  lower_volc_thickness * lower_volc_CO2 + \
                  transition_thickness * transition_CO2 + \
                  sheeted_dykes_thickness * sheeted_dykes_CO2 + \
                  gabbros_thickness * CO2_gabbro)/(thickness * volcanic_percent/100)

              if CO2_volcanic < 0:
                  CO2_volcanic == 0

          max_DS = 100
          m = carbon_max/max_DS
          #print tmp_DS
          if velocity > 40:
               #for fast ridges just assume max C is .32? questionable, no data available i think
               carbon_max = .32
          m = carbon_max/max_DS
          #check potential mistake in this line

          CO2_serpentinites_vertical_area = thickness * \
                                 asymmetry_factor * \
                                 per * 1/100 * \
                                 1000.0/865.0 * \
                                 tmp_DS * m * 1/100 *   velocity

          CO2_volcanics_vertical_area = thickness * \
                               asymmetry_factor * \
                               volcanic_percent * 1/100 * \
                               CO2_volcanic * 1/100 * velocity#* \

          CO2_thickness_serpentinites = thickness * \
                                 per * 1/100 * \
                                 1000.0/865.0 * \
                                 tmp_DS * m *1/100

          CO2_thickness_volcanics = thickness * \
                               volcanic_percent * 1/100 * \
                               CO2_volcanic * 1/100#* \

          serpentinites_vertical_area = thickness * \
                                 asymmetry_factor * \
                                 per * 1/100 * \
                                 1000.0/865.0 * \
                                 tmp_DS * m * velocity

          serpentinites_point_thickness = thickness * \
                                 per * 1/100 * \
                                 1000.0/865.0 * \
                                 tmp_DS * m

          peridotites_vertica_area = thickness * \
                                 asymmetry_factor * \
                                 per * 1/100 * velocity

          peridotites_point_thickness = thickness * \
                                 per * 1/100

          gabbros_vertica_area = gabbros_thickness * \
                                 asymmetry_factor * velocity

          sheeted_dykes_vertical_area = sheeted_dykes_thickness * \
                                 asymmetry_factor * velocity

          transition_vertical_area = transition_thickness * \
                                 asymmetry_factor * velocity

          lower_volcs_vertical_area = lower_volc_thickness * \
                                 asymmetry_factor * velocity

          upper_volcs_vertical_area = upper_volc_thickness * \
                                 asymmetry_factor * velocity


          CO2_thickness = CO2_thickness_serpentinites + CO2_thickness_volcanics
          CO2_vertical_area = CO2_serpentinites_vertical_area + CO2_volcanics_vertical_area

        #print DS,upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2,CO2_gabbro
          C_storage[time][SR].append(velocity)
          C_storage[time][per].append(per)
          C_storage[time][DS].append((tmp_DS,upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2,CO2_gabbro))
          C_storage[time][thick].append((thickness,upper_volc_thickness, lower_volc_thickness,transition_thickness,sheeted_dykes_thickness,gabbros_thickness, peridotites_point_thickness))
          C_storage[time][vertical_area].append((upper_volcs_vertical_area, lower_volcs_vertical_area,transition_vertical_area,sheeted_dykes_vertical_area,gabbros_vertica_area, peridotites_vertica_area))
          C_storage[time][C_total].append(CO2_thickness*12./44.)#multiply by mass C/mass CO2 to get C
          C_storage[time][C_total_vertical].append(CO2_vertical_area*12./44.)
          C_storage[time][C_volcanic].append((CO2_thickness_volcanics, CO2_volcanic)) #wt%
          C_storage[time][C_serpentinites].append((CO2_thickness_serpentinites, (tmp_DS * m)))
          C_storage[time][serpentinites_thick].append(serpentinites_point_thickness)
          C_storage[time][serpentinites_vertical].append(serpentinites_vertical_area)
          C_storage[time][Index].append(index)
          C_storage[time][length].append(row.Length_km)
          C_storage[time][boundary].append(row.FeatureType)
          C_storage[time][StartLat].append(row.StartPointLat)
          C_storage[time][StartLon].append(row.StartPointLon)
          C_storage[time][EndLat].append(row.EndPointLat)
          C_storage[time][EndLon].append(row.EndPointLon)
          C_storage[time][RightPlateID].append(row.RightPlate)
          C_storage[time][LeftPlateID].append(row.LeftPlate)


  #save if wanted

  #date = datetime.today().strftime('%Y-%m-%d')
  #filename = 'C_storage_%s_pacific_test.p' % date
  #outfile = open('%s/%s' % (loaddir, filename), 'wb')
  #pickle.dump(C_storage, outfile)
  #outfile.close()

  return(filename)
