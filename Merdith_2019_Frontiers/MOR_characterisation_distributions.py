import numpy as np
from scipy.stats import truncnorm

def spreading_rate_generator():

    '''
    returns a random spreading rate between 0 and 100 mm/a defined by the
    preserved seafloor spreading rate record (e.g. Muller et al. 2016)
    '''

     #possible range of spreading rate
    SR_range = (1,101,1)
    SR_mu = 37 #mean spreading rate after Muller et al. 2016
    SR_sigma = 27 #standard deviation of spreading rate after Muller et al., 2016
    SR_length = len(SR_range)
    #this line creates a normal distribution between two values (SR max and SR min)
    #with a shifted mean and standard deviation
    SR_normal = truncnorm((min(SR_range)-SR_mu)/SR_sigma,
                              (max(SR_range)-SR_mu)/SR_sigma,
                                 loc=SR_mu, scale=SR_sigma)
    spreading_rate = SR_normal.rvs() #generate random number

    return spreading_rate

def SR_and_thickness_slow_ultraslow(spreading_rate):
    '''
    this function takes a spreading rate and returns a thickness for slow or ultraslow spreading
    (less than 40 mm/a)
    and returns a 'thickness' which in this case is the maximum depth of water penetration giving the
    maximum depth of serpentinisation.
    i.e. this returns the depth to the unserpentinised mantle lithosphere
    '''
    #print spreading_rate
    if 0 <= spreading_rate <= 20:
        #Ultraslow thickness proportion normal distribution (to unaltered mantle lithosphere)
        thickness_range = np.arange(3,4.2,.2)
    elif 20 < spreading_rate <= 40:
        #Slow thickness proportion normal distribution (to unaltered mantle lithosphere)
        thickness_range = np.arange(3,4.2,.2)

    thickness_mu = np.mean(thickness_range)
    thickness_sigma = np.std(thickness_range)
    thickness_length = len(thickness_range)
    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
                              (max(thickness_range)-thickness_mu)/thickness_sigma,
                             loc=thickness_mu, scale=thickness_sigma)
    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
    thickness = np.random.choice(thickness_normal)

    return thickness

def SR_and_thickness_inter_fast(DS, volcanic_percent):
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

    upper_volc_thickness = np.random.choice(np.arange(upper_volc[1]-upper_volc[1]*0.2, upper_volc[1]+upper_volc[1]*0.2 + upper_volc[1]*0.04, upper_volc[1]*0.04))
    lower_volc_thickness = np.random.choice(np.arange(lower_volc[1]-lower_volc[1]*0.2, lower_volc[1]+lower_volc[1]*0.2 + lower_volc[1]*0.04, lower_volc[1]*0.04))
    transition_thickness = np.random.choice(np.arange(transition[1]-transition[1]*0.2, transition[1]+transition[1]*0.2 + transition[1]*0.04, transition[1]*0.04))
    sheeted_dykes_thickness = np.random.choice(np.arange(sheeted_dykes[1]-sheeted_dykes[1]*0.2, sheeted_dykes[1]+sheeted_dykes[1]*0.2 + sheeted_dykes[1]*0.04, sheeted_dykes[1]*0.04))
    gabbros_thickness = np.random.choice(np.arange(gabbros[1]-gabbros[1]*0.2, gabbros[1]+gabbros[1]*0.2 + gabbros[1]*0.04, gabbros[1]*0.04))

    total_volc_thickness = upper_volc_thickness + lower_volc_thickness + transition_thickness + sheeted_dykes_thickness + gabbros_thickness

    total_thickness = total_volc_thickness/volcanic_percent * 100

    peridotite_thickness_tmp = total_thickness - total_volc_thickness

    peridotite_thickness = peridotite_thickness_tmp - (peridotite_thickness_tmp * DS/100)

    serpentinite_thickness = peridotite_thickness_tmp * DS/100 #volume change?* 1000/865)**(1./3/) #get the added z factor of thickness to serpentites

    new_total_thickness = total_volc_thickness + peridotite_thickness + serpentinite_thickness

    return new_total_thickness, upper_volc_thickness, lower_volc_thickness, transition_thickness, sheeted_dykes_thickness, gabbros_thickness

def degree_of_serpentinisation(dsSurf, thick):
    '''
    this function takes the degree of serpentinisation at the surface of the
    oceanic lithosphere column and returns the total degree of serpentinisation
    for the rock body only for lithosphere exhumed at slow and ultraslow ridges
    '''
    if dsSurf == 0:
        serp_total = 0
    else:
        range_top = np.arange(.8,1.45,.05)
        mu, sigma, length = np.mean(range_top), np.std(range_top), len(range_top)
        dTop_normal = truncnorm((min(range_top)-mu)/sigma,
                                  (max(range_top)-mu)/sigma,
                                 loc=mu, scale=sigma)
        dTop_normal = dTop_normal.rvs(length) #generate random numbers
        dTop = np.random.choice(dTop_normal)
        dBot = thick #as below is unaltered mantle peridotites
        #print 'check', dBot, dTop, dsSurf
        area_total = dBot*100 #(max area of dsSurf)
        serp_area = (dsSurf * dTop) + (((dBot - dTop) * dsSurf)/2)

        serp_total = serp_area/area_total * 100

    return serp_total

def SR_and_dsSurf_slow_ultraslow(spreading_rate, thickness):

    thickness=thickness
    possible_dsSurf = []

    SR =int(spreading_rate) #also slip rate for transform faults

    #initiate some SR
    if SR <= 20:
        X =  np.linspace(1,20,21)
        y = np.arange(100,79,-1)
        y1 = 0.0449*SR**2 - 1.899*SR + 98.47
    if 20 < SR <= 40:
        X = np.linspace(21,40,81)
        y = np.arange(100,19,-1)
        y1 = 0.2193*SR**2 - 17.63*SR + 369.03

    ind, val = min(enumerate(X), key=lambda x: abs(x[1]-SR))
    if y1 < 0:
        y1 = 0

    possible_dsSurf.append(np.arange(np.round(y1), y[ind]+1,1))

    for j in possible_dsSurf:
        mu, sigma, length = np.mean(j), np.std(j), len(j)

        tmp = truncnorm((min(j)-mu)/sigma,
                       (max(j)-mu)/sigma,
                       loc=mu, scale=sigma)
        tmp = tmp.rvs(length)
        dsSurf = np.random.choice(tmp)

    DS = degree_of_serpentinisation(dsSurf, thickness)

    return DS

def SR_and_dsSurf_inter_fast(spreading_rate):

    possible_dsSurf = []

    SR =int(spreading_rate) #also slip rate for transform faults
    if SR > 100:
        SR = 100

    #initiate some SR
    if SR <= 70:
        X = np.linspace(41,70,21)
        y = np.arange(20,-1,-1)
        y1 = 0.0231*SR**2 - 3.247*SR + 112.77
    if 70 < SR:
        X = np.linspace(71,100,11)
        y = np.arange(10,-1,-1)
        y1 = 0.0115*SR**2 - 2.3162*SR + 115.48

    ind, val = min(enumerate(X), key=lambda x: abs(x[1]-SR))
    if y1 < 0:
        y1 = 0

    possible_dsSurf.append(np.arange(np.round(y1), y[ind]+2,1)) #+2 to avoid an index error at the extreme values

    for j in possible_dsSurf:
        mu, sigma, length = np.mean(j), np.std(j), len(j)

        tmp = truncnorm((min(j)-mu)/sigma,
                       (max(j)-mu)/sigma,
                       loc=mu, scale=sigma)
        tmp = tmp.rvs(length)
        dsSurf = np.random.choice(tmp)

    DS = dsSurf

    return DS

def SR_and_thickness(spreading_rate):
    #print spreading_rate
    if 0 <= spreading_rate <= 20:
        #Ultraslow thickness proportion normal distribution (to unaltered mantle lithosphere)
        thickness_range = np.arange(3,4.2,.2)
    elif 20 < spreading_rate <= 40:
        #Slow thickness proportion normal distribution (to unaltered mantle lithosphere)
        thickness_range = np.arange(3,4.2,.2)
    elif 40 < spreading_rate <= 70:
        #Intermediate thickness proportion normal distribution
        thickness_range = np.arange(6,7.2,.2)
    elif 70 < spreading_rate:
        #Fast thickness proportion normal distribution
        thickness_range = np.arange(6,7.2,.2)

    thickness_mu = np.mean(thickness_range)
    thickness_sigma = np.std(thickness_range)
    thickness_length = len(thickness_range)
    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
                              (max(thickness_range)-thickness_mu)/thickness_sigma,
                             loc=thickness_mu, scale=thickness_sigma)
    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
    thickness = np.random.choice(thickness_normal)
    return thickness

def SR_and_peridotite(spreading_rate):

    if 0 <= spreading_rate <= 20:
        #Ultraslow peridotite proportion normal distribution
        peridotite_range = np.arange(80,101,1)
    elif 20 < spreading_rate <= 40:
        #Slow peridotite proportion normal distribution
        peridotite_range = np.arange(12.5,80,1)
    elif 40 < spreading_rate <= 70:
        #Intermediate peridotite proportion normal distribution
        peridotite_range = np.arange(5,15,1)
    elif 70 < spreading_rate:
        #Fast peridotite proportion normal distribution
        peridotite_range = np.arange(0,10,1)


    peridotite_mu = np.mean(peridotite_range)
    peridotite_sigma = np.std(peridotite_range)
    peridotite_length = len(peridotite_range)
    peridotite_normal = truncnorm((min(peridotite_range)-peridotite_mu)/peridotite_sigma,
                              (max(peridotite_range)-peridotite_mu)/peridotite_sigma,
                             loc=peridotite_mu, scale=peridotite_sigma)
    peridotite_normal = peridotite_normal.rvs(peridotite_length) #generate random numbers
    peridotite = np.random.choice(peridotite_normal)

    return peridotite

def SR_and_transforms(spreading_rate):
    peridotite_range = np.arange(0,101,1)
    thickness_range = np.arange(6,7.2,.2)

    peridotite_mu = np.mean(peridotite_range)
    peridotite_sigma = np.std(peridotite_range)
    peridotite_length = len(peridotite_range)
    peridotite_normal = truncnorm((min(peridotite_range)-peridotite_mu)/peridotite_sigma,
                              (max(peridotite_range)-peridotite_mu)/peridotite_sigma,
                             loc=peridotite_mu, scale=peridotite_sigma)
    peridotite_normal = peridotite_normal.rvs(peridotite_length) #generate random numbers
    peridotite = np.random.choice(peridotite_normal)

    thickness_mu = np.mean(thickness_range)
    thickness_sigma = np.std(thickness_range)
    thickness_length = len(thickness_range)
    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
                              (max(thickness_range)-thickness_mu)/thickness_sigma,
                             loc=thickness_mu, scale=thickness_sigma)
    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
    thickness = np.random.choice(thickness_normal)

    return peridotite, thickness

def carbon_content_slow_ultraslow(spreading_rate, peridotite_percent):

    #figure out how much CO2 in volcanics preserved in slow ridges (i.e. 100% - peridotite)

    CO2_gabbro_mu = 0.096
    CO2_gabbro_sigma = 0.05
    CO2_gabbro = float(np.random.normal(CO2_gabbro_mu, CO2_gabbro_sigma, 1))
    if CO2_gabbro < 0:
        CO2_gabbro = 0
    #carbon initial from Kelemen and Manning, the amount of carbon in 100% +f serp (decreasing to 0% at 0 DS)
    carbon_initial_range = np.arange(0.32,0.37,0.01)
    carbon_initial = np.random.choice(carbon_initial_range)

    return carbon_initial, CO2_gabbro

def carbon_content_inter_fast(spreading_rate, peridotite_percent, bottom_water_temp):

    if 5 > bottom_water_temp >= 0:
        bottom_water_temperature_multiplier = 0.4
    elif 10 > bottom_water_temp >= 5:
        bottom_water_temperature_multiplier = 0.6
    elif 15 > bottom_water_temp >= 10:
        bottom_water_temperature_multiplier = 0.8
    elif 20 > bottom_water_temp >= 15:
        bottom_water_temperature_multiplier = 1
    elif bottom_water_temp >= 20:
        bottom_water_temperature_multiplier = 1.2


    #bottom water controls only for upper 300 m of oceanic crust(i.e. upper volcs and lower volcs)
    #Gills and Coogan, 2011

    upper_volc = 2.5 * bottom_water_temperature_multiplier
    lower_volc = 0.171 * bottom_water_temperature_multiplier
    transition = 0.073
    sheeted_dykes = 0.145
    gabbros = 0.096

    upper_volc_CO2 = np.random.choice(np.arange(upper_volc-upper_volc*0.2, upper_volc+upper_volc*0.2 + upper_volc*0.04, upper_volc*0.04))
    lower_volc_CO2 = np.random.choice(np.arange(lower_volc-lower_volc*0.2, lower_volc+lower_volc*0.2 + lower_volc*0.04, lower_volc*0.04))
    transition_CO2 = np.random.choice(np.arange(transition-transition*0.2, transition+transition*0.2 + transition*0.04, transition*0.04))
    sheeted_dykes_CO2 = np.random.choice(np.arange(sheeted_dykes-sheeted_dykes*0.2, sheeted_dykes+sheeted_dykes*0.2 + sheeted_dykes*0.04, sheeted_dykes*0.04))

    CO2_gabbro_mu = 0.096
    CO2_gabbro_sigma = 0.05
    CO2_gabbro = float(np.random.normal(CO2_gabbro_mu, CO2_gabbro_sigma, 1))

    return upper_volc_CO2, lower_volc_CO2, transition_CO2, sheeted_dykes_CO2, CO2_gabbro, bottom_water_temperature_multiplier
