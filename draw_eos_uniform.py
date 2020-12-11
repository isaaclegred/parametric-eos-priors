# Test Sampling eos's from parametric values
import numpy as np
import scipy.interpolate as interp
import lalsimulation as lalsim
import lal
import argparse
from matplotlib import pyplot as plt
import astropy.constants as const
import astropy.units as u


c = const.c.cgs.value
# Characteristic refinement number (probably should be changed on a case by case basis)
N = 100

# Model of equation of state prior to sample from: 
# Need to sample gamma1, gamma2, gamma3, and logp1
# From the Lackey and Wade paper (I don't actually )
logp1_range = (33.531, 34.5)
gamma1_range  = (1.52, 5)
gamma2_range = (1.05, 5)
gamma3_range = (1, 5)

parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws")

# need
class eos_polytrope:
    def __init__(self,logp1, gamma1, gamma2, gamma3):
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3
        self.logp1 = logp1
        self.eos = lalsim.SimNeutronStarEOS4ParameterPiecewisePolytrope(
                                                             self.logp1-1,
                                                            self.gamma1, 
                                                            self.gamma2, 
                                                            self.gamma3)
        print(logp1, gamma1, gamma2, gamma3)
        self.family = lalsim.CreateSimNeutronStarFamily(self.eos)
    # Get the eos family from the paramaters. 
    def get_eos(self):
        return self.family.eos
    # Evaluate the eos in terms of epsilon(p)
    def eval_energy_density(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):    
            eps = np.zeros(len(p))  
            for i, pres in enumerate(p):
                eps[i] = lalsim.SimNeutronStarEOSEnergyDensityOfPressure(pres,self.eos)    
        else:
             eps = lalsim.SimNeutronStarEOSEnergyDensityOfPressure(p, self.eos)
        return eps
    def eval_phi(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):    
            eps = np.zeros(len(p))  
            for i, pres in enumerate(p):
                eps[i] = lalsim.SimNeutronStarEOSEnergyDensityDerivOfPressure(pres,self.eos)    
        else:
             eps = lalsim.SimNeutronStarEOSEnergyDensityDerivOfPressure(p, self.eos)
        return eps
    def eval_baryon_density(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):    
            rho = np.zeros(len(p))  
            for i, pres in enumerate(p):
                rho [i] = lalsim.SimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(
                    lalsim.SimNeutronStarEOSPseudoEnthalpyOfPressure(pres,self.eos), self.eos)    
        else:
            rho  = lalsim.SimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(
                lalsim.SimNeutronStarEOSPseudoEnthalpyOfPressure(p,self.eos), self.eos) 
        return rho
    
# We also need gamma1 > gamma2 > gamma3 for thermodynamic stability, so I guess we sample gamma3 first 
# and then constrain the others based on this. This is the (somewhat) uniform prior on the gamma's, I think
# I still need to glue together the 
def get_eos_realization_uniform_poly (logp1_range, gamma1_range, gamma2_range, gamma3_range):
    eps = .1
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(gamma2_range[0]+eps, gamma1 - eps)
    gamma3 = np.random.uniform(gamma3_range[0], gamma2 - eps) 
    logp1 = np.random.uniform(*logp1_range) 
    return eos_polytrope(logp1, gamma1, gamma2, gamma3)    
 
# Because we have an analytic equation of state, we can compute the derivative dmu/dp 
# analytically.  Therefore we can compute phi analytically (Doesn't seem to actually be necessary)


# Stitch EoS onto the known EoS below nuclear saturation density. 
# Use Sly log(p1) = 34.384, gamma1 = 3.005, gamma2 = 2.988, gamma3 = 2.851
# There's some subtlety here related to where the pressure is known.  Here
# it is given at the dividing line between rho_1 and rho_2,  but we want it
# at the stitching point (I think this is fine though, because we can just 
# evaluate the modelt at the stiching point)

#SLy model
sly_polytrope_model = eos_polytrope(34.384, 3.005, 2.988, 2.851)

# It's probably the case that from this point on Universality can handle it, (hopefully)...

 
def create_eos_draw_file(name):
    eos_poly = get_eos_realization_uniform_poly(logp1_range, gamma1_range, gamma2_range, gamma3_range)
    # FIXME WORRY ABOUT CGS VS SI!!!!! (Everything is in SI till the last step :/ ) 
    p_sly = np.geomspace(1.0e31, 3.9e32, 100)
    p_main = np.geomspace (3.9e33, 9.0e37, 800)
    eps_sly = sly_polytrope_model.eval_energy_density(p_sly)
    eps_main = eos_poly.eval_energy_density(p_main)
    rho_b_sly = eos_poly.eval_baryon_density(p_sly)
    rho_b_main = eos_poly.eval_baryon_density(p_main)
    p = np.concatenate([p_sly, p_main])
    eps = np.concatenate([eps_sly, eps_main])
    rho_b = np.concatenate([rho_b_sly, rho_b_main])
    data = np.transpose(np.stack([p/c**2*10 , eps/c**2*10, rho_b])) # *10 because Everything above is done in SI
    np.savetxt(name,data, header = 'pressurec2, energy_densityc2, baryon_density', fmt='%.10e')

if __name__ == "__main__":
    args = parser.parse_args()
    num_draws = args.num_draws
    for i in range(num_draws):
        name = "eos-draw-" + "%06d" %i + ".csv"
        create_eos_draw_file(name)
