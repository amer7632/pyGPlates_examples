import numpy as np
import pandas as pd


import MOR_characterisation_serp_flux_FINAL


#select random pacific spreading rate
def random_spreading(samples,distribution):
    y = np.random.randint(1001,size=samples)
    return distribution[y]

def crust_characterisation(start, stop, step, df, samples, bottom_water_temperature_curve):
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
    prop_per = 'Peridotite'
    DS = 'Degree of Serpentinisation'
    CO2wt_mean = 'CO2 mean'
    CO2wt_std = 'CO2 std'
    thick_mean = 'Thickness Mean'
    thick_std = 'Thickness std'
    vertical_area_mean = 'Vertical Area Mean'
    vertical_area_std = 'Vertical Area std'
    C_total = 'C-prod'
    C_total_vertical = 'C-prod_vertical'
    C_volcanic_serp = 'C-Volcanic_serp'
    serpentinites_thick = 'Serpentinites_raw'
    serpentinites_vertical = 'Serpentinites_raw_vertical'
    PDF_params = 'PDF Parameters'
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
        C_storage[time][prop_per] = []
        C_storage[time][DS] = []
        C_storage[time][CO2wt_mean] = []
        C_storage[time][CO2wt_std] = []
        C_storage[time][thick_mean] = []
        C_storage[time][thick_std] = []
        C_storage[time][vertical_area_mean] = []
        C_storage[time][vertical_area_std] = []
        C_storage[time][serpentinites_thick] = []
        C_storage[time][serpentinites_vertical] = []
        C_storage[time][C_volcanic_serp] = []
        C_storage[time][C_total] = []
        C_storage[time][C_total_vertical] = []
        C_storage[time][PDF_params] = []
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
            velocity = np.ones(samples)
            velocity = velocity*row.Plate_Velocity*10

           #block out if not using plate model
            if time > 180:
                new_conjugate_plates = (row.RightPlate,row.LeftPlate)
                if all(x in panthalassa_plate_ids for x in new_conjugate_plates) == True:
                    velocity = np.ones(samples)
                    velocity = velocity*random_spreading(sample,spread_dist)*10
                    #print 'here', velocity

            #calculate the mean bottom water temperature for the first 20 Ma existence of a parcel of ocean crust
            temp_tmp = []
            bottom_water_temperature = np.ones(samples)
            for i in np.arange(time,time-21,-1):
                if i < 0:
                    break
                temp_tmp.append(bottom_water_temperature_curve[i])
            bottom_water_temperature = bottom_water_temperature * np.mean(temp_tmp)

            per = MOR_characterisation_serp_flux_FINAL.SR_and_peridotite(samples,
                                                                         velocity)

            volcanic_percent = 100 - per
            asymmetry_factor = 1
            calc_length = row.Length_km

            #calculate variables for mid ocean ridge segments

            calc_length = row.Length_km
            width = velocity
            #as all velocitie are the same for each point, we just check the first
            if velocity[0] <= 40:
                asymmetry_factor = 1
                thickness = MOR_characterisation_serp_flux_FINAL.SR_and_thickness_slow_ultraslow(samples)

                tmp_DS = MOR_characterisation_serp_flux_FINAL.SR_and_dsSurf_slow_ultraslow(samples,
                                                                                           velocity,
                                                                                          thickness) #DS and thickness is called within this function
                #print DS
                carbon_max, CO2_gabbro = MOR_characterisation_serp_flux_FINAL.carbon_content_slow_ultraslow(samples,
                                                                                                            velocity)
                upper_volc_CO2 = lower_volc_CO2 = transition_CO2 = sheeted_dykes_CO2 = 0
                upper_volc_thickness = lower_volc_thickness = transition_thickness = sheeted_dykes_thickness = 0
                gabbros_thickness = thickness * volcanic_percent/100

                #incase a neg CO2 value comes back (i think this was fixed in previous tweaks, but
                #leaving it in to be sure)
                indices_neg_CO2_gabbro = CO2_gabbro < 0
                CO2_gabbro[indices_neg_CO2_gabbro] = 0
                CO2_volcanic = gabbros_thickness * CO2_gabbro

            else:
                carbon_max = 0
                tmp_DS = MOR_characterisation_serp_flux_FINAL.SR_and_dsSurf_inter_fast(samples,
                                                                                   velocity) #DS and thickness is called within this function

                thickness, upper_volc_thickness, lower_volc_thickness, transition_thickness, \
                sheeted_dykes_thickness, gabbros_thickness \
                = MOR_characterisation_serp_flux_FINAL.SR_and_thickness_inter_fast(samples,tmp_DS,volcanic_percent)

                upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2, \
                CO2_gabbro, bottom_water_temperature_multiplier = MOR_characterisation_serp_flux_FINAL.carbon_content_inter_fast(samples,
                                                                                            velocity,
                                                                                            bottom_water_temperature)
                #incase a neg CO2 value comes back (i think this was fixed in previous tweaks, but
                #leaving it in to be sure)
                indices_neg_CO2_gabbro = CO2_gabbro < 0
                CO2_gabbro[indices_neg_CO2_gabbro] = 0

                #gives us total CO2 storage in volcanics as a fraction (NOT A PERCENT)
                CO2_volcanic =  (upper_volc_thickness * upper_volc_CO2 + \
                    lower_volc_thickness * lower_volc_CO2 + \
                    transition_thickness * transition_CO2 + \
                    sheeted_dykes_thickness * sheeted_dykes_CO2 + \
                    gabbros_thickness * CO2_gabbro)/(thickness * volcanic_percent/100)

                #incase a neg CO2 value comes back (i think this was fixed in previous tweaks, but
                #leaving it in to be sure)
                indices_neg_CO2_volcanic = CO2_volcanic < 0
                CO2_volcanic[indices_neg_CO2_volcanic] = 0


            #carbon in serpentinite
            max_DS = 100
            #print tmp_DS
            if velocity[0] > 40:
                 #for fast ridges just assume max C is .32? questionable, no data available i think
                carbon_max = .32
            m = carbon_max/max_DS

            C_thickness_volcanics = thickness * volcanic_percent * 1/100 * CO2_volcanic * 1/100 * 12./44.
            C_thickness_serpentinites = thickness * per * 1/100 * 1000.0/865.0 * tmp_DS * m *1/100 * 12./44.
            C_total_volc_serp =  C_thickness_volcanics + C_thickness_serpentinites

            peridotites_point_thickness = thickness * per * 1/100
            serpentinites_point_thickness = peridotites_point_thickness * 1000.0/865.0 * tmp_DS *1/100

            #print PDF_parameters
            C_storage[time][SR].append(np.mean(velocity))
            C_storage[time][prop_per].append(np.mean(per))
            C_storage[time][DS].append(np.mean(tmp_DS))
            C_storage[time][CO2wt_mean].append((np.mean(upper_volc_CO2),
                                                 np.mean(lower_volc_CO2),
                                                 np.mean(transition_CO2),
                                                 np.mean(sheeted_dykes_CO2),
                                                 np.mean(CO2_gabbro)))
            C_storage[time][thick_mean].append((np.mean(thickness),
                                                np.mean(upper_volc_thickness),
                                                np.mean(lower_volc_thickness),
                                                np.mean(transition_thickness),
                                                np.mean(sheeted_dykes_thickness),
                                                np.mean(gabbros_thickness),
                                                np.mean(peridotites_point_thickness)))
            C_storage[time][vertical_area_mean].append((np.mean(upper_volc_thickness * velocity),
                                                        np.mean(lower_volc_thickness * velocity),
                                                        np.mean(transition_thickness * velocity),
                                                        np.mean(sheeted_dykes_thickness * velocity),
                                                        np.mean(gabbros_thickness * velocity),
                                                        np.mean(peridotites_point_thickness * velocity)))
            #multiply by mass C/mass CO2 to get C
            C_storage[time][C_total].append(np.mean(C_total_volc_serp))

            C_storage[time][C_total_vertical].append(np.mean(C_total_volc_serp)*velocity)

            C_storage[time][C_volcanic_serp].append((np.mean(C_thickness_volcanics),
                                                     np.mean(CO2_volcanic * 12./44.),
                                                     np.mean(C_thickness_serpentinites),
                                                     np.mean(tmp_DS * m)))

            C_storage[time][serpentinites_thick].append(np.mean(serpentinites_point_thickness))

            C_storage[time][serpentinites_vertical].append(np.mean(serpentinites_point_thickness) * velocity)

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

    date = datetime.today().strftime('%Y-%m-%d')
    filename = 'C_storage_%s_plate_model_CHAPMAN.p' % date
    outfile = open('%s/%s' % (loaddir, filename), 'wb')
    pickle.dump(C_storage, outfile)
    outfile.close()

    return(filename)
