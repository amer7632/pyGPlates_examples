import numpy as np
from scipy.stats import truncnorm

def SR_and_dsSurf_slow_ultraslow(samples, spreading_rate, thickness, normal=True):

	'''
	Return degree of serpentinisation given a spreading rate. Optional flag is
	for normal or uniform distribution. Default is normal. Only for slow and
	ultraslow spreading ridges (<= 40 mm/a)
	'''

	# with vectors, all spreading rates are the same

	spreading_rate_round = np.round(spreading_rate[0])
	if spreading_rate_round == 0:
		DS_slow_ultraslow = np.zeros(samples)

	else:
		if spreading_rate_round <= 20:
			X = np.linspace(1,20,21) #if SR <= 20:
			y = np.arange(100,79,-1) # min depending on spreading rate
			y1 = 0.0449*spreading_rate_round**2 - 1.899*spreading_rate_round + 98.47 #max
			if y1<0:
				y1=0
			ind = (np.round(spreading_rate_round/0.95)-1).astype(int)
		else:
			X = np.linspace(21,40,81) #if 20 < SR <= 40:
			y = np.arange(100,19,-1)
			y1 = 0.2193*spreading_rate_round**2 - 17.63*spreading_rate_round + 369.03
			if y1<0:
					y1=0

			ind = (np.round((spreading_rate_round-21)/0.2375)-1).astype(int)
			if ind < 0:
				ind = 0

		yi = ind
		if normal:
		#print 'normal'
			mu = np.mean([np.round(y1),y[yi]])
			sigma = np.std(np.arange(np.round(y1), y[yi]+2,1))
			dsSurf = truncnorm.rvs((np.round(y1)-mu)/sigma,
			   (y[yi]-mu)/sigma,
			   loc=mu, scale=sigma,size=samples)

			DS_slow_ultraslow = DS_and_thickness(samples, dsSurf, thickness, True)

		else:
			#print 'uniform'

			dsSurf = np.random.uniform(y1, y[yi], size=samples)

			DS_slow_ultraslow = DS_and_thickness(samples, dsSurf, thickness, False)

	return DS_slow_ultraslow

def SR_and_dsSurf_inter_fast(samples, spreading_rate, normal=True):

	'''
Return degree of serpentinisation given a spreading rate. Optional flag is
for normal or uniform distribution. Default is normal. only for intermediate and
fast spreading ridges (> 40 mm/a)
	'''

	# with vectors for conditional probability
	#possible_dsSurf = []
	spreading_rate_round = np.round(spreading_rate[0])

	if spreading_rate_round > 100: spreading_rate_round = 100
	if spreading_rate_round == 40: spreading_rate_round = 41

	if spreading_rate_round <=70:
		y = np.arange(20,-1,-1) # min DS depending on spreading rate
		y1 = 0.0231*spreading_rate_round**2 - 3.247*spreading_rate_round + 112.77
		if y1 < 0:
			y1 = 0
		ind = np.round((spreading_rate_round-41)/1.45)-1
		#to make sure always a 0 or above index
		if ind < 0:
			ind = 0
	else:
		y = np.arange(10,-1,-1)
		y1 = 0.0115*spreading_rate_round**2 - 2.3162*spreading_rate_round + 115.48
		if y1 < 0:
			y1 = 0
		ind = np.round((spreading_rate_round-71)/2.9)-1
		if ind < 0:
			ind = 0

	yi = int(ind)
	if normal:
		#print 'normal'
		mu = np.mean([np.round(y1),y[yi]+1])
		sigma = np.std(np.arange(np.round(y1), y[yi]+2,1))
		DS_inter_fast = truncnorm.rvs((np.round(y1)-mu)/sigma,
		   (y[yi]-mu)/sigma,
		   loc=mu, scale=sigma,size=samples)

	else:
		#print 'uniform'
		DS_inter_fast = np.random.uniform(y1, y[yi], size=samples)

	return DS_inter_fast

def DS_and_thickness(samples, dsSurf, thickness, normal=True):
	'''
	This function takes a degree of serpentinisation and thickness
	and returns the total degree of serpentinisation for the system.
	Optional flag is for normal or uniform distribution. Default is normal.
	'''

	range_top = np.arange(.8,1.45,.05)
	if normal:
		#print 'normal'
		mu, sigma, length = np.mean(range_top), np.std(range_top), len(range_top)
		dInflex_tmp = truncnorm((min(range_top)-mu)/sigma,
		  (max(range_top)-mu)/sigma,
		 loc=mu, scale=sigma)
		dInflex_tmp = dInflex_tmp.rvs(length) #generate random numbers
		dInflex = np.random.choice(dInflex_tmp, size=samples)

	else:
		#print 'uniform'
		dInflex = np.random.uniform(0.8, 1.4, size=samples)


	dBot = thickness #as below is unaltered mantle peridotites
	area_total = dBot*100 #(max area of dsSurf)

	serp_area = (np.abs(dsSurf) * dInflex) + (((dBot - dInflex) * np.abs(dsSurf))/2.)

	serp_total = serp_area/area_total * 100

	return serp_total

def SR_and_thickness(samples, spreading_rate, normal=True):

	'''
Return thickness given a spreading rate. Optional flag is for normal or uniform
distribution. Default is normal. This is just for bulk crustal thicknesses, without
subdivisions (as in the H2 paper)
	'''

	# with vectors
	thickness = np.zeros(samples)

	SR1 = np.where(spreading_rate <= 40)
	SR2 = np.where(spreading_rate > 40)

	#truncated normal distributions for thickness
	if normal:
		thickness = truncnorm.rvs(-3.75,3.75, scale=0.2,size=samples) # mean=0, std=0.4, min/max=3.75/0.4=1.5
		thickness[SR1] = thickness[SR1]+3.5
		thickness[SR2] = thickness[SR2]+6.5

	#uniform distributions
	else:
		#print 'uniform'
		if len(SR1[0]) != 0:
			thickness[SR1] = np.random.uniform(3.0, 4.0, size=samples)
		if len(SR2[0]) != 0:
			thickness[SR2] = np.random.uniform(6.0, 8.0, size=samples)

	return thickness

def SR_and_thickness_slow_ultraslow(samples):
    # with vectors
    '''
    this function takes a spreading rate and returns a thickness for slow or
	ultraslow spreading (less than 40 mm/a) and returns a 'thickness' which in
	this case is the maximum depth of water penetration giving the maximum
	depth of serpentinisation. i.e. this returns the depth to the
	 unserpentinised mantle lithosphere
    '''
    #truncated normal distributions for thickness
    thickness = truncnorm.rvs(-0.6,0.6, scale=0.4,size=samples) +3.6

    return thickness

def SR_and_thickness_inter_fast(samples, DS, volcanic_percent):
	# with vectors
	'''
	this function takes a spreading rate and other tectonic parameters of oceanic crust produced
	and returns the thickness of the different volcanic layers (basalts, dykes, gabbros etc.)
	'''
	#proportions of normal lithosphere
	#[carbon%, thickness, thickness% of total depth[]
	total_depth = 7.0
	upper_volc = [2.5, 0.3, 0.3/total_depth]
	lower_volc = [0.171, 0.3, 0.3/total_depth]
	transition = [0.073, 0.2, 0.2/total_depth]
	sheeted_dykes = [0.145, 1.2, 1.2/total_depth]
	gabbros = [0.096, 5, 5/total_depth] #stdev 0.05

	upper_volc_thickness = np.random.choice(np.arange(upper_volc[1]-upper_volc[1]*0.2, upper_volc[1]+upper_volc[1]*0.2 + upper_volc[1]*0.04, upper_volc[1]*0.04),size=samples)
	lower_volc_thickness = np.random.choice(np.arange(lower_volc[1]-lower_volc[1]*0.2, lower_volc[1]+lower_volc[1]*0.2 + lower_volc[1]*0.04, lower_volc[1]*0.04),size=samples)
	transition_thickness = np.random.choice(np.arange(transition[1]-transition[1]*0.2, transition[1]+transition[1]*0.2 + transition[1]*0.04, transition[1]*0.04),size=samples)
	sheeted_dykes_thickness = np.random.choice(np.arange(sheeted_dykes[1]-sheeted_dykes[1]*0.2, sheeted_dykes[1]+sheeted_dykes[1]*0.2 + sheeted_dykes[1]*0.04, sheeted_dykes[1]*0.04),size=samples)
	gabbros_thickness = np.random.choice(np.arange(gabbros[1]-gabbros[1]*0.2, gabbros[1]+gabbros[1]*0.2 + gabbros[1]*0.04, gabbros[1]*0.04), size=samples)

	total_volc_thickness = upper_volc_thickness + lower_volc_thickness + transition_thickness + sheeted_dykes_thickness + gabbros_thickness

	total_thickness = total_volc_thickness/volcanic_percent * 100

	#peridotite thickness before any degree_of_serpentinisation
	peridotite_thickness_tmp = total_thickness - total_volc_thickness

	#after serpentinisation
	peridotite_thickness = peridotite_thickness_tmp - (peridotite_thickness_tmp * DS/100.)

	serpentinite_thickness = peridotite_thickness_tmp * DS/100. #volume change?* 1000/865)**(1./3/) #get the added z factor of thickness to serpentites

	new_total_thickness = total_volc_thickness + peridotite_thickness + serpentinite_thickness

	return new_total_thickness, upper_volc_thickness, lower_volc_thickness, transition_thickness, sheeted_dykes_thickness, gabbros_thickness

def SR_and_peridotite(samples, spreading_rate, normal=True):

	'''
	Return peridotite proportion given a spreading rate. Optional flag is for
	normal or uniform distribution. Default is normal.
	'''

	# with vectors
	#print len(spreading_rate)
	peridotite = np.zeros(samples)

	SR1 = np.where(spreading_rate<=20)
	SR2 = np.where((spreading_rate>20) & (spreading_rate<=40))
	SR3 = np.where((spreading_rate>40) & (spreading_rate<=70))
	SR4 = np.where(spreading_rate>70)

	if normal:
	    #print 'normal'
	    if len(SR1[0]) != 0:
	        peridotite[SR1] = truncnorm.rvs(-1.65,1.65, scale=6,size=samples) + 90
	    if len(SR2[0]) != 0:
	        peridotite[SR2] = truncnorm.rvs(-1.71,1.71, scale=19.6,size=samples) + 46
	    if len(SR3[0]) != 0:
	        peridotite[SR3] = truncnorm.rvs(-1.57,1.57, scale=2.87,size=samples) + 9.5
	    if len(SR4[0]) != 0:
	        peridotite[SR4] = truncnorm.rvs(-1.57,1.57, scale=2.87,size=samples) + 5

	else:
	    #print 'uniform'
	    if len(SR1[0]) != 0:
	        peridotite[SR1] = np.random.uniform(80, 100, size=samples)
	    if len(SR2[0]) != 0:
	        peridotite[SR2] = np.random.uniform(12.5, 80, size=samples)
	    if len(SR3[0]) != 0:
	        peridotite[SR3] = np.random.uniform(5, 15, size=samples)
	    if len(SR4[0]) != 0:
	        peridotite[SR4] = np.random.uniform(0, 10, size=samples)

	return peridotite

def SR_and_transforms(samples, spreading_rate, normal=True):

	'''
Return thickness and peridotite given a slip rate. Optional flag is for normal or uniform
distribution. Default is normal.
	'''

	# with vectors

	if normal:
		#print 'normal'
		peridotite = truncnorm.rvs(-1.71,1.71, scale=29.2,size=samples) +50
		thickness = truncnorm.rvs(-1.5,1.5, scale=0.4,size=samples) +6.6
		DS = truncnorm.rvs(-2.5,1.5, scale=20,size=samples) + 50

	else:
		#print 'uniform'
		peridotite = np.random.uniform(0, 100, size=samples)
		thickness = np.random.uniform(6, 8, size=samples)
		DS = np.random.uniform(0, 80, size=samples)

	return peridotite, thickness, DS

def carbon_content_slow_ultraslow(samples, spreading_rate):

    '''
    figure out how much CO2 in volcanics preserved in slow ridges (i.e. 100% - peridotite)
    we just use gabbro for the moment, the distribution given in Kelemen and Manning (2015) is normal
    so no option for uniform.
    '''
    # with vectors

    CO2_gabbro_mu = 0.096
    CO2_gabbro_sigma = 0.05
    CO2_gabbro = np.random.normal(CO2_gabbro_mu, CO2_gabbro_sigma, size=samples)
    CO2_gabbro[CO2_gabbro<0] = 0

    #carbon initial from Kelemen and Manning, the amount of carbon in 100% +f serp (decreasing to 0% at 0 DS)
    carbon_initial_range = np.arange(0.32,0.37,0.01)
    carbon_initial = np.random.choice(carbon_initial_range,size=samples) ## SEA_random

    return carbon_initial, CO2_gabbro

def carbon_content_inter_fast(samples, spreading_rate, bottom_water_temp):
    # with vectors

    T1 = np.where(bottom_water_temp<5)
    T2 = np.where((bottom_water_temp>=5) & (bottom_water_temp<10))
    T3 = np.where((bottom_water_temp>=10) & (bottom_water_temp<15))
    T4 = np.where((bottom_water_temp>=15) & (bottom_water_temp<20))
    T5 = np.where(bottom_water_temp>=20)

    bottom_water_temperature_multiplier = np.ones(samples)
    bottom_water_temperature_multiplier[T1] = 0.4
    bottom_water_temperature_multiplier[T2] = 0.6
    bottom_water_temperature_multiplier[T3] = 0.8
    bottom_water_temperature_multiplier[T4] = 1
    bottom_water_temperature_multiplier[T5] = 1.2

    #bottom water controls only for upper 300 m of oceanic crust(i.e. upper volcs and lower volcs)
    #Gills and Coogan, 2011
	#these are in wt%
    upper_volc = 2.5 * bottom_water_temperature_multiplier
    lower_volc = 0.171 * bottom_water_temperature_multiplier
    transition = 0.073
    sheeted_dykes = 0.145
    gabbros = 0.096

    upper_volc_CO2 = (np.random.rand(samples)+(upper_volc-0.5))*(upper_volc*0.4)
    lower_volc_CO2 = (np.random.rand(samples)+(lower_volc-0.5))*(lower_volc*0.4)
    transition_CO2 = (np.random.rand(samples)+(transition-0.5))*(transition*0.4)
    sheeted_dykes_CO2 = (np.random.rand(samples)+(sheeted_dykes-0.5))*(sheeted_dykes*0.4)

    CO2_gabbro_mu = 0.096
    CO2_gabbro_sigma = 0.05
    CO2_gabbro = np.random.normal(CO2_gabbro_mu, CO2_gabbro_sigma, size=samples)

    return upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2, CO2_gabbro, bottom_water_temperature_multiplier

def carbon_content_serp(samples, carbon_max, spreading_rate):

    max_DS = 100.

    C1 = np.where(spreading_rate<=40)
    C2 = np.where(spreading_rate>40)

    CO2_serp = np.zeros(samples)

    CO2_serp[C1] = carbon_max/max_DS
    CO2_serp[C2] = 0.32

    return CO2_serp

def get_random_Fe3_wt(degree_of_serpentinisation, samples):
    '''
    '''
    #print degree_of_serpentinisation, samples
    curveFit = np.asarray([0., 0.07452397, 0.14851618, 0.22197664, 0.29490535, 0.36730231, 0.43916752, 0.51050098, 0.58130268, 0.65157263,
       0.72131083, 0.79051728, 0.85919197, 0.92733492, 0.99494611, 1.06202555, 1.12857324, 1.19458918, 1.26007336, 1.32502579,
       1.38944648, 1.4533354 , 1.51669258, 1.57951801, 1.64181168, 1.7035736 , 1.76480377, 1.82550219, 1.88566886, 1.94530377,
       2.00440694, 2.06297835, 2.12101801, 2.17852591, 2.23550207, 2.29194647, 2.34785912, 2.40324002, 2.45808917, 2.51240657,
       2.56619221, 2.6194461 , 2.67216825, 2.72435863, 2.77601727, 2.82714416, 2.87773929, 2.92780267, 2.9773343 , 3.02633418,
       3.0748023 , 3.12273868, 3.1701433 , 3.21701617, 3.26335729, 3.30916665, 3.35444427, 3.39919013, 3.44340424, 3.4870866 ,
       3.53023721, 3.57285607, 3.61494317, 3.65649852, 3.69752212, 3.73801397, 3.77797406, 3.81740241, 3.856299  , 3.89466384,
       3.93249693, 3.96979827, 4.00656785, 4.04280569, 4.07851177, 4.1136861 , 4.14832868, 4.1824395 , 4.21601858, 4.2490659 ,
       4.28158147, 4.31356529, 4.34501735, 4.37593767, 4.40632623, 4.43618304, 4.4655081 , 4.49430141, 4.52256296, 4.55029277,
       4.57749082, 4.60415712, 4.63029167, 4.65589447, 4.68096551, 4.7055048 , 4.72951234, 4.75298813, 4.77593217, 4.79834445,
       4.82022499])
    fitError =np.asarray([0., 0.0126386 , 0.02528133, 0.03793232, 0.0505957 , 0.06327558, 0.07597606, 0.08870126, 0.10145524, 0.11424209,
       0.12706584, 0.13993054, 0.15284019, 0.16579878, 0.17881027, 0.19187859, 0.20500765, 0.21820131, 0.23146343, 0.2447978 ,
       0.25820819, 0.27169834, 0.28527194, 0.29893264, 0.31268406, 0.32652977, 0.34047329, 0.35451811, 0.36866765, 0.38292532,
       0.39729444, 0.41177833, 0.42638021, 0.44110329, 0.45595071, 0.47092558, 0.48603093, 0.50126975, 0.51664501, 0.53215958,
       0.54781631, 0.56361799, 0.57956736, 0.5956671 , 0.61191985, 0.6283282 , 0.64489469, 0.66162179, 0.67851194, 0.69556752,
       0.71279088, 0.7301843 , 0.74775001, 0.76549021, 0.78340704, 0.80150259, 0.81977892, 0.83823802, 0.85688187, 0.87571236,
       0.89473138, 0.91394075, 0.93334225, 0.95293763, 0.9727286, 0.9927168 , 1.01290387, 1.0332914 , 1.05388092, 1.07467394,
       1.09567193, 1.11687634, 1.13828856, 1.15990995, 1.18174186, 1.20378557, 1.22604235, 1.24851344, 1.27120004, 1.29410333,
       1.31722445, 1.34056451, 1.36412461, 1.38790579, 1.4119091, 1.43613553, 1.46058607, 1.48526167, 1.51016327, 1.53529177,
       1.56064805, 1.58623298, 1.61204738, 1.63809209, 1.6643679, 1.69087557, 1.71761587, 1.74458954, 1.77179728, 1.79923981,
       1.82691779])

    rounded_DS = np.round(degree_of_serpentinisation, 0)
    rounded_DS = rounded_DS.astype(int)

    upper_bound = np.zeros([samples])
    lower_bound = np.zeros([samples])
    fe3_range = np.zeros([samples])
    fe3_values = np.zeros([samples])

    upper_bound += curveFit[rounded_DS] + fitError[rounded_DS]
    lower_bound += curveFit[rounded_DS] - fitError[rounded_DS]

    fe3_range = np.array([np.linspace(i,j,10) for i,j in zip(lower_bound,upper_bound)])

    for ind, i in enumerate(fe3_range):

        fe3_values[ind] = np.random.choice(fe3_range[ind])

    return fe3_values

def ferric_iron_moles_proportion_per_ds(degree_of_serpentinisation, samples):
	'''
	This function returns the molar prorportion of ferric iron in serpentine
	in serpentinite after Fig. 8c in Andreani et al. 2013 (Lithos)
	'''
	ferric_iron_mole_proportion = np.asarray([0.0, 0.52250336, 0.58682945, 0.63336639, 0.66848787,
				0.69584329, 0.71767263, 0.73542678, 0.75008673, 0.76233946,
				0.77268028, 0.78147519, 0.78900044, 0.79546831, 0.80104458,
				0.80586046, 0.81002114, 0.81361179, 0.81670211, 0.81934955,
				0.82160188, 0.82349904, 0.82507462, 0.82635705, 0.82737042,
				0.82813529, 0.8286692 , 0.82898716, 0.82910206, 0.82902489,
				0.8287651 , 0.82833075, 0.82772868, 0.82696471, 0.8260437 ,
				0.82496972, 0.82374606, 0.82237538, 0.82085974, 0.81920063,
				0.81739908, 0.81545563, 0.81337042, 0.81114317, 0.80877327,
				0.80625972, 0.80360122, 0.80079615, 0.79784258, 0.79473832,
				0.79148088, 0.78806751, 0.78449519, 0.78076066, 0.77686041,
				0.77279068, 0.76854747, 0.76412655, 0.75952346, 0.75473352,
				0.74975179, 0.74457316, 0.73919226, 0.73360354, 0.72780122,
				0.72177933, 0.71553169, 0.70905194, 0.70233354, 0.69536977,
				0.68815374, 0.68067844, 0.67293668, 0.66492118, 0.65662455,
				0.64803932, 0.63915793, 0.62997283, 0.62047642, 0.61066117,
				0.60051958, 0.59004425, 0.57922795, 0.5680636 , 0.55654442,
				0.54466388, 0.53241583, 0.51979458, 0.5067949 , 0.49341218,
				0.47964243, 0.46548247, 0.45092993, 0.4359834 , 0.42064252,
				0.40490808, 0.38878217, 0.37226824, 0.35537124, 0.33809777,
				0.32045614])

	DS = np.round(degree_of_serpentinisation, 0)
	DS = DS.astype(int)
	ferric_iron_mole_propotion_in_serpentine = np.zeros([samples])

    #print DS
	ferric_iron_mole_propotion_in_serpentine += ferric_iron_mole_proportion[DS]
	return ferric_iron_mole_propotion_in_serpentine
