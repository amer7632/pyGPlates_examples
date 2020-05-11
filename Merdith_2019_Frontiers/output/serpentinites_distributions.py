import numpy as np
from scipy.stats import truncnorm


def serpentinite_distribution_vertical(samples):
    #set up distributions
    mu_vertical = 0.05343248058082458 #np.mean(serpentinites_vertical) taken from PDF notebook
    sigma_vertical =  0.033612691919305444#np.std(serpentinites_vertical) taken from PDF notebook
    max_vertical = 0.17067791805475954 #max(serpentinites_vertical) taken from PDF notebook
    
    vertical_range = np.linspace(0,max_vertical,50) #possible range 

    #this line creates a normal distribution between two values (serpentinite max and serpentinite min)
    vertical_normal = truncnorm((min(vertical_range)-mu_vertical)/sigma_vertical, 
                              (max(vertical_range)-mu_vertical)/sigma_vertical,
                                 loc=mu_vertical, scale=sigma_vertical)
    
    vertical_perturbation = vertical_normal.rvs(size=samples)
    
    vertical_ratio = vertical_perturbation/mu_vertical
    return vertical_ratio
   
def serpentinite_distribution_point(samples):#set up distributions
   
    mu_point = 0.002533120763009994 #np.mean(serpentinites_vertical) taken from PDF notebook
    sigma_point = 0.0026442364129314107 #np.std(serpentinites_vertical) taken from PDF notebook
    max_point = 0.010924536325850663 #max(serpentinites_vertical) taken from PDF notebook

    point_range = np.linspace(0.,max_point,50) #possible range 
   
    #this line creates a normal distribution between two values (serpentinite max and serpentinite min)
    #with a shifted mean and standard deviation
    point_normal = truncnorm((min(point_range)-mu_point)/sigma_point, 
                              (max(point_range)-mu_point)/sigma_point,
                                 loc=mu_point, scale=sigma_point)
    
    point_perturbation = point_normal.rvs(size=samples)
    point_ratio = point_perturbation/mu_point
    return point_ratio