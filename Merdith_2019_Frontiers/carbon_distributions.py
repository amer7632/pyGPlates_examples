import numpy as np
from scipy.stats import truncnorm

def carbon_distribution_vertical(samples):
    #set up distributions
    mu_vertical = 0.738886948862488 #np.mean(CO2_vertical) taken from PDF notebook
    sigma_vertical =  0.8164297218821015 #np.std(CO2_vertical) taken from PDF notebook
    max_vertical = 4.671950891575566 #max(CO2_vertical) taken from PDF notebook
    
    vertical_range = np.linspace(0,max_vertical,50) #possible range 

    #this line creates a normal distribution between two values (CO2 max and CO2 min)
    vertical_normal = truncnorm((min(vertical_range)-mu_vertical)/sigma_vertical, 
                              (max(vertical_range)-mu_vertical)/sigma_vertical,
                                 loc=mu_vertical, scale=sigma_vertical)
    
    vertical_perturbation = vertical_normal.rvs(size=samples)
    
    vertical_ratio = vertical_perturbation/mu_vertical
    return vertical_ratio
   
def carbon_distribution_point(samples):#set up distributions
   
    mu_point = 0.013260060933599575 #np.mean(CO2_final) taken from PDF notebook
    sigma_point = 0.012878892516125031 #np.std(CO2_final) taken from PDF notebook
    max_point = 00.06300102532689246 #max(CO2_final) taken from PDF notebook

    point_range = np.linspace(0.,max_point,50) #possible range 
   
    #this line creates a normal distribution between two values (CO2 max and CO2 min)
    #with a shifted mean and standard deviation
    point_normal = truncnorm((min(point_range)-mu_point)/sigma_point, 
                              (max(point_range)-mu_point)/sigma_point,
                                 loc=mu_point, scale=sigma_point)
    
    point_perturbation = point_normal.rvs(size=samples)
    point_ratio = point_perturbation/mu_point
    return point_ratio


