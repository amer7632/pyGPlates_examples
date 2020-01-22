import numpy as np
from scipy.stats import truncnorm


def SR_and_thickness_slow_ultraslow(spreading_rate):
    # with vectors
    '''
    this function takes a spreading rate and returns a thickness for slow or ultraslow spreading
    (less than 40 mm/a)
    and returns a 'thickness' which in this case is the maximum depth of water penetration giving the
    maximum depth of serpentinisation.
    i.e. this returns the depth to the unserpentinised mantle lithosphere
    '''
##    #print spreading_rate
##    if 0 <= spreading_rate <= 20:
##        #Ultraslow thickness proportion normal distribution (to unaltered mantle lithosphere)
##        thickness_range = np.arange(3,4.2,.2)
##    elif 20 < spreading_rate <= 40:
##        #Slow thickness proportion normal distribution (to unaltered mantle lithosphere)
##        thickness_range = np.arange(3,4.2,.2)
##
##    thickness_mu = np.mean(thickness_range)
##    thickness_sigma = np.std(thickness_range)
##    thickness_length = len(thickness_range)
##    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
##                              (max(thickness_range)-thickness_mu)/thickness_sigma,
##                             loc=thickness_mu, scale=thickness_sigma)
##    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
##    thickness = np.random.choice(thickness_normal) ## SEA_random


    #truncated normal distributions for thickness
    peridotite = np.zeros([len(spreading_rate)])
    thickness = truncnorm.rvs(-0.6,0.6, scale=0.4,size=len(spreading_rate)) +3.6

    return thickness

def SR_and_thickness_inter_fast(DS, volcanic_percent):
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

    upper_volc_thickness = np.random.choice(np.arange(upper_volc[1]-upper_volc[1]*0.2, upper_volc[1]+upper_volc[1]*0.2 + upper_volc[1]*0.04, upper_volc[1]*0.04),len(DS))
    lower_volc_thickness = np.random.choice(np.arange(lower_volc[1]-lower_volc[1]*0.2, lower_volc[1]+lower_volc[1]*0.2 + lower_volc[1]*0.04, lower_volc[1]*0.04),len(DS))
    transition_thickness = np.random.choice(np.arange(transition[1]-transition[1]*0.2, transition[1]+transition[1]*0.2 + transition[1]*0.04, transition[1]*0.04),len(DS))
    sheeted_dykes_thickness = np.random.choice(np.arange(sheeted_dykes[1]-sheeted_dykes[1]*0.2, sheeted_dykes[1]+sheeted_dykes[1]*0.2 + sheeted_dykes[1]*0.04, sheeted_dykes[1]*0.04),len(DS))
    gabbros_thickness = np.random.choice(np.arange(gabbros[1]-gabbros[1]*0.2, gabbros[1]+gabbros[1]*0.2 + gabbros[1]*0.04, gabbros[1]*0.04))

    total_volc_thickness = upper_volc_thickness + lower_volc_thickness + transition_thickness + sheeted_dykes_thickness + gabbros_thickness

    total_thickness = total_volc_thickness/volcanic_percent * 100

    peridotite_thickness_tmp = total_thickness - total_volc_thickness

    peridotite_thickness = peridotite_thickness_tmp - (peridotite_thickness_tmp * DS/100)

    serpentinite_thickness = peridotite_thickness_tmp * DS/100 #volume change?* 1000/865)**(1./3/) #get the added z factor of thickness to serpentites

    new_total_thickness = total_volc_thickness + peridotite_thickness + serpentinite_thickness

    return new_total_thickness, upper_volc_thickness, lower_volc_thickness, transition_thickness, sheeted_dykes_thickness, gabbros_thickness

def DS_and_thickness(dsSurf, thick):
    # with vectors
    '''
    this function takes the degree of serpentinisation at the surface of the oceanic lithosphere column
    and returns the total degree of serpentinisation for the rock body
    only for lithosphere exhumed at slow and ultraslow ridges
    '''
    serp_total = np.zeros([len(thick)])

##    range_top = np.arange(.8,1.45,.05)
##    mu, sigma, length = np.mean(range_top), np.std(range_top), len(range_top)
##    dTop_normal = truncnorm((min(range_top)-mu)/sigma,
##                              (max(range_top)-mu)/sigma,
##                             loc=mu, scale=sigma)
##    dTop_normal = dTop_normal.rvs(length) #generate random numbers
##    dTop = np.random.choice(dTop_normal) ## SEA_random
    dTop = truncnorm.rvs(-1.6,1.6,scale=0.19,size=len(thick))+1.1
    dBot = thick #as below is unaltered mantle peridotites
    #print 'check', dBot, dTop, dsSurf
    area_total = dBot*100 #(max area of dsSurf)
    serp_area = (dsSurf * dTop) + (((dBot - dTop) * dsSurf)/2)

    serp_total = serp_area/area_total * 100
    serp_total[dsSurf==0]=0

    return serp_total

def SR_and_dsSurf_slow_ultraslow(spreading_rate, thickness):
    # with vectors

    possible_dsSurf = []

    #also slip rate for transform faults
    spreading_rate[spreading_rate>100] = 100

    SR1 = np.where(spreading_rate<=20) # assumes can't get more than 40 in
    SR2 = np.where(spreading_rate>20)

    ind = np.zeros([len(spreading_rate)])
    y1 = np.zeros([len(spreading_rate)])

    XS = np.linspace(1,20,21) #if SR <= 20:
    yS = np.arange(100,79,-1) # min depending on spreading rate
    y1[SR1] = 0.0449*spreading_rate[SR1]**2 - 1.899*spreading_rate[SR1] + 98.47 #max

    ind[SR1] = np.round(spreading_rate[SR1]/0.95)-1

    XF = np.linspace(21,40,81) #if 20 < SR <= 40:
    yF = np.arange(100,19,-1)
    y1[SR2] = 0.2193*spreading_rate[SR2]**2 - 17.63*spreading_rate[SR2] + 369.03
    y1[y1<0]=0

    ind[SR2] = (np.round((spreading_rate[SR2]-21)/0.2375)-1)
    ind[ind<0] = 0

 #   ind, val = min(enumerate(X), key=lambda x: abs(x[1]-SR)) # returns the index and value of hte closest point in X to SR


 #   possible_dsSurf.append(np.arange(np.round(y1), y[ind]+1,1))
    dsSurf=np.zeros(len(spreading_rate))
    DS=np.zeros(len(spreading_rate))

    for i in range(0,len(SR1[0])):
        j=SR1[0][i]
        yi = ind[j].astype(int)
        mu = np.mean([np.round(y1[j]),yS[yi]])
        sigma = np.std(np.arange(np.round(y1[j]), yS[yi]+1,1))
        dsSurf[j] = truncnorm.rvs((np.round(y1[j])-mu)/sigma,
                       (yS[yi]-mu)/sigma,
                       loc=mu, scale=sigma)

    for i in range(0,len(SR2[0])):
        j=SR2[0][i]
        yi = ind[j].astype(int)
        mu = np.mean([np.round(y1[j]),yF[yi]])
        sigma = np.std(np.arange(np.round(y1[j]), yF[yi]+1,1))
        dsSurf[j] = truncnorm.rvs((np.round(y1[j])-mu)/sigma,
                       (yF[yi]-mu)/sigma,
                       loc=mu, scale=sigma)

    DS = DS_and_thickness(dsSurf, thickness)


    return DS

def SR_and_dsSurf_slow_ultraslow_CP(spreading_rate, thickness):
    # with vectors, all spreading rates are the same

    possible_dsSurf = []

    #also slip rate for transform faults

    if spreading_rate[0]<=20:
        X = np.linspace(1,20,21) #if SR <= 20:
        y = np.arange(100,79,-1) # min depending on spreading rate
        y1 = 0.0449*spreading_rate[0]**2 - 1.899*spreading_rate[0] + 98.47 #max
        if y1<0:
            y1=0
        ind = (np.round(spreading_rate[0]/0.95)-1).astype(int)
    else:
        X = np.linspace(21,40,81) #if 20 < SR <= 40:
        y = np.arange(100,19,-1)
        y1 = 0.2193*spreading_rate[0]**2 - 17.63*spreading_rate[0] + 369.03
        if y1<0:
            y1=0

        ind = (np.round((spreading_rate[0]-21)/0.2375)-1).astype(int)
        if ind<0:
            ind = 0



 #   ind, val = min(enumerate(X), key=lambda x: abs(x[1]-SR)) # returns the index and value of hte closest point in X to SR


 #   possible_dsSurf.append(np.arange(np.round(y1), y[ind]+1,1)

    yi = ind
    mu = np.mean([np.round(y1),y[yi]])
    sigma = np.std(np.arange(np.round(y1), y[yi]+1,1))
    dsSurf = truncnorm.rvs((np.round(y1)-mu)/sigma,
                   (y[yi]-mu)/sigma,
                   loc=mu, scale=sigma,size=len(spreading_rate))


    DS = DS_and_thickness(dsSurf, thickness)


    return DS

def SR_and_dsSurf_inter_fast(spreading_rate):
    # with vectors

    possible_dsSurf = []

##    SR =int(spreading_rate) #also slip rate for transform faults
    spreading_rate[spreading_rate>100] = 100

##    #initiate some SR
##    if SR <= 70:
##        X = np.linspace(41,70,21)
##        y = np.arange(20,-1,-1)
##        y1 = 0.0231*SR**2 - 3.247*SR + 112.77
##    if 70 > SR:
##        X = np.linspace(71,100,11)
##        y = np.arange(10,-1,-1)
##        y1 = 0.0115*SR**2 - 2.3162*SR + 115.48
##
##    ind, val = min(enumerate(X), key=lambda x: abs(x[1]-SR))
##    if y1 < 0:
##        y1 = 0

    SR1 = np.where(spreading_rate<=70) # assumes can't get more than 40 in
    SR2 = np.where(spreading_rate>70)

    ind = np.zeros([len(spreading_rate)])
    y1 = np.zeros([len(spreading_rate)])

    yS = np.arange(20,-1,-1) # min depending on spreading rate
    y1[SR1] = 0.0231*spreading_rate[SR1]**2 - 3.247*spreading_rate[SR1] + 112.77 #max

    ind[SR1] = np.round((spreading_rate[SR1]-41)/1.45)-1

    yF = np.arange(10,-1,-1)
    y1[SR2] = 0.0115*spreading_rate[SR2]**2 - 2.3162*spreading_rate[SR2] + 115.48
    y1[y1<0]=0

    ind[SR2] = (np.round((spreading_rate[SR2]-71)/2.9)-1).astype(int)
    ind[ind<0] = 0

    #possible_dsSurf.append(np.arange(np.round(y1), y[ind]+2,1)) #+2 to avoid an index error at the extreme values

    DS=np.zeros(len(spreading_rate))

    for i in range(0,len(SR1[0])):
        j=SR1[0][i]
        yi = ind[j].astype(int)
        mu = np.mean([np.round(y1[j]),yS[yi]+1])
        sigma = np.std(np.arange(np.round(y1[j]), yS[yi]+2,1))
        DS[j] = truncnorm.rvs((np.round(y1[j])-mu)/sigma,
                       (yS[yi]-mu)/sigma,
                       loc=mu, scale=sigma)
    for i in range(0,len(SR2[0])):
        j=SR2[0][i]
        yi = ind[j].astype(int)
        mu = np.mean([np.round(y1[j]),yF[yi]+1])
        sigma = np.std(np.arange(np.round(y1[j]), yF[yi]+2,1))
        DS[j] = truncnorm.rvs((np.round(y1[j])-mu)/sigma,
                       (yF[yi]-mu)/sigma,
                       loc=mu, scale=sigma)

##    for j in possible_dsSurf:
##        mu, sigma, length = np.mean(j), np.std(j), len(j)
##
##        tmp = truncnorm((min(j)-mu)/sigma,
##                       (max(j)-mu)/sigma,
##                       loc=mu, scale=sigma)
##        tmp = tmp.rvs(length)
##        dsSurf = np.random.choice(tmp)  ## SEA_random

##    DS = dsSurf

    return DS

def SR_and_dsSurf_inter_fast_CP(spreading_rate):
    # with vectors for conditional probability

    possible_dsSurf = []

##    SR =int(spreading_rate) #also slip rate for transform faults
    spreading_rate[spreading_rate>100] = 100

    SR1 = np.where(spreading_rate<=70) # assumes can't get more than 40 in
    SR2 = np.where(spreading_rate>70)

    if spreading_rate[0] <=70:
        y = np.arange(20,-1,-1) # min depending on spreading rate
        y1 = 0.0231*spreading_rate[0]**2 - 3.247*spreading_rate[0] + 112.77 #max
        ind = np.round((spreading_rate[0]-41)/1.45)-1
        if ind < 0:
            ind = 0
    else:
        y = np.arange(10,-1,-1)
        y1 = 0.0115*spreading_rate[0]**2 - 2.3162*spreading_rate[0] + 115.48
        if y1 < 0:
            y1=0
        ind = np.round((spreading_rate[0]-71)/2.9)-1
        if ind < 0:
            ind = 0

    yi = int(ind)
    mu = np.mean([np.round(y1),y[yi]+1])
    sigma = np.std(np.arange(np.round(y1), y[yi]+2,1))
    DS = truncnorm.rvs((np.round(y1)-mu)/sigma,
                   (y[yi]-mu)/sigma,
                   loc=mu, scale=sigma,size=len(spreading_rate))

    return DS


def SR_and_thickness(spreading_rate):
    # with vectors

    SR1 = np.where(spreading_rate<=40)
    SR2 = np.where(spreading_rate>40)

    #truncated normal distributions for thickness
    thickness = truncnorm.rvs(-3.75,3.75, scale=0.4,size=len(spreading_rate)) # mean=0, std=0.4, min/max=3.75/0.4=1.5
    thickness[SR1] = thickness[SR1]+3.6
    thickness[SR2] = thickness[SR2]+6.6

##
##    #print spreading_rate
##    if 0 <= spreading_rate <= 20:
##        #Ultraslow thickness proportion normal distribution (to unaltered mantle lithosphere)
##        thickness_range = np.arange(3,4.2,.2)
##    elif 20 < spreading_rate <= 40:
##        #Slow thickness proportion normal distribution (to unaltered mantle lithosphere)
##        thickness_range = np.arange(3,4.2,.2)
##    elif 40 < spreading_rate <= 70:
##        #Intermediate thickness proportion normal distribution
##        thickness_range = np.arange(6,7.2,.2)
##    elif 70 < spreading_rate:
##        #Fast thickness proportion normal distribution
##        thickness_range = np.arange(6,7.2,.2)
##
##    thickness_mu = np.mean(thickness_range)
##    thickness_sigma = np.std(thickness_range)
##    thickness_length = len(thickness_range)
##    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
##                              (max(thickness_range)-thickness_mu)/thickness_sigma,
##                             loc=thickness_mu, scale=thickness_sigma)
##    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
##    thickness = np.random.choice(thickness_normal)

    return thickness

def SR_and_peridotite(spreading_rate):
    # with vectors

    SR1 = np.where(spreading_rate<=20)
    SR2 = np.where((spreading_rate>20) & (spreading_rate<=40))
    SR3 = np.where((spreading_rate>40) & (spreading_rate<=70))
    SR4 = np.where(spreading_rate>70)

    #truncated normal distributions for thickness
    peridotite = np.zeros([len(spreading_rate)])
    peridotite[SR1] = truncnorm.rvs(-1.65,1.65, scale=6,size=len(SR1[0])) +90
    peridotite[SR2] = truncnorm.rvs(-1.71,1.71, scale=19.6,size=len(SR2[0])) +46
    peridotite[SR3] = truncnorm.rvs(-1.57,1.57, scale=2.87,size=len(SR3[0])) +9.5
    peridotite[SR4] = truncnorm.rvs(-1.57,1.57, scale=2.87,size=len(SR4[0])) +5

##    if 0 <= spreading_rate <= 20:
##        #Ultraslow peridotite proportion normal distribution
##        peridotite_range = np.arange(80,101,1)
##    elif 20 < spreading_rate <= 40:
##        #Slow peridotite proportion normal distribution
##        peridotite_range = np.arange(12.5,80,1)
##    elif 40 < spreading_rate <= 70:
##        #Intermediate peridotite proportion normal distribution
##        peridotite_range = np.arange(5,15,1)
##    elif 70 < spreading_rate:
##        #Fast peridotite proportion normal distribution
##        peridotite_range = np.arange(0,10,1)
##
##
##    peridotite_mu = np.mean(peridotite_range)
##    peridotite_sigma = np.std(peridotite_range)
##    peridotite_length = len(peridotite_range)
##    peridotite_normal = truncnorm((min(peridotite_range)-peridotite_mu)/peridotite_sigma,
##                              (max(peridotite_range)-peridotite_mu)/peridotite_sigma,
##                             loc=peridotite_mu, scale=peridotite_sigma)
##    peridotite_normal = peridotite_normal.rvs(peridotite_length) #generate random numbers
##    peridotite = np.random.choice(peridotite_normal) ## SEA_random

    return peridotite

def SR_and_transforms(spreading_rate):
    # with vectors

##    peridotite_range = np.arange(0,101,1)
##    thickness_range = np.arange(6,7.2,.2)
##
##    peridotite_mu = np.mean(peridotite_range)
##    peridotite_sigma = np.std(peridotite_range)
##    peridotite_length = len(peridotite_range)
##    peridotite_normal = truncnorm((min(peridotite_range)-peridotite_mu)/peridotite_sigma,
##                              (max(peridotite_range)-peridotite_mu)/peridotite_sigma,
##                             loc=peridotite_mu, scale=peridotite_sigma)
##    peridotite_normal = peridotite_normal.rvs(peridotite_length) #generate random numbers
##    peridotite = np.random.choice(peridotite_normal) ## SEA_random
    peridotite = truncnorm.rvs(-1.71,1.71, scale=29.2,size=len(spreading_rate)) +50

##    thickness_mu = np.mean(thickness_range)
##    thickness_sigma = np.std(thickness_range)
##    thickness_length = len(thickness_range)
##    thickness_normal = truncnorm((min(thickness_range)-thickness_mu)/thickness_sigma,
##                              (max(thickness_range)-thickness_mu)/thickness_sigma,
##                             loc=thickness_mu, scale=thickness_sigma)
##    thickness_normal = thickness_normal.rvs(thickness_length) #generate random numbers
##    thickness = np.random.choice(thickness_normal) ## SEA_random
    thickness = truncnorm.rvs(-1.5,1.5, scale=0.4,size=len(spreading_rate)) +6.6

    return peridotite, thickness

def carbon_content_slow_ultraslow(spreading_rate, peridotite_percent):
    # with vectors

    #figure out how much CO2 in volcanics preserved in slow ridges (i.e. 100% - peridotite)

    CO2_gabbro_mu = 0.096
    CO2_gabbro_sigma = 0.05
    CO2_gabbro = np.random.normal(CO2_gabbro_mu, CO2_gabbro_sigma, len(spreading_rate))
    CO2_gabbro[CO2_gabbro<0] = 0

    #carbon initial from Kelemen and Manning, the amount of carbon in 100% +f
    #serp (decreasing to 0% at 0 DS)
    carbon_initial_range = np.arange(0.32,0.37,0.01)
    carbon_initial = np.random.choice(carbon_initial_range,len(spreading_rate)) ## SEA_random

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
