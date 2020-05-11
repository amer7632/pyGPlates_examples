import numpy as np
from scipy.stats import truncnorm


def volcanic_uncertainty_thickness_fast(samples):

    #we assumed 20% uncertainty on thickness layers

    fast_volcanic_perturbation = np.random.uniform(0.8, 1.2,size=samples)

    return fast_volcanic_perturbation

def serps_uncertainty_thickness_fast(samples):

    mu_vertical = 75 #np.mean(serpentinites_vertical) taken from PDF notebook
    sigma_vertical =  61#np.std(serpentinites_vertical) taken from PDF notebook
    max_vertical = 310 #max(serpentinites_vertical) taken from PDF notebook

    vertical_range = np.linspace(0,max_vertical,50) #possible range

    #this line creates a normal distribution between two values (serpentinite max and serpentinite min)
    vertical_normal = truncnorm((min(vertical_range)-mu_vertical)/sigma_vertical,
                              (max(vertical_range)-mu_vertical)/sigma_vertical,
                                 loc=mu_vertical, scale=sigma_vertical)

    vertical_perturbation = vertical_normal.rvs(size=samples)

    fast_serp_perturbation = vertical_perturbation/mu_vertical
    return fast_serp_perturbation

def volcanic_uncertainty_thickness_slow(samples):

    mu_vertical = 1318 #np.mean(serpentinites_vertical) taken from PDF notebook
    sigma_vertical =  907#np.std(serpentinites_vertical) taken from PDF notebook
    max_vertical = 3612 #max(serpentinites_vertical) taken from PDF notebook

    vertical_range = np.linspace(0,max_vertical,50) #possible range

    #this line creates a normal distribution between two values (serpentinite max and serpentinite min)
    vertical_normal = truncnorm((min(vertical_range)-mu_vertical)/sigma_vertical,
                              (max(vertical_range)-mu_vertical)/sigma_vertical,
                                 loc=mu_vertical, scale=sigma_vertical)

    vertical_perturbation = vertical_normal.rvs(size=samples)

    slow_volc_perturbation = vertical_perturbation/mu_vertical
    return slow_volc_perturbation

def serps_uncertainty_thickness_slow(samples):

    mu_vertical = 1236 #np.mean(serpentinites_vertical) taken from PDF notebook
    sigma_vertical = 799#np.std(serpentinites_vertical) taken from PDF notebook
    max_vertical = 2989 #max(serpentinites_vertical) taken from PDF notebook

    vertical_range = np.linspace(0,max_vertical,50) #possible range

    #this line creates a normal distribution between two values (serpentinite max and serpentinite min)
    vertical_normal = truncnorm((min(vertical_range)-mu_vertical)/sigma_vertical,
                              (max(vertical_range)-mu_vertical)/sigma_vertical,
                                 loc=mu_vertical, scale=sigma_vertical)

    vertical_perturbation = vertical_normal.rvs(size=samples)

    slow_serp_perturbation = vertical_perturbation/mu_vertical
    return slow_serp_perturbation
