#Test Sampling eos's from spectral coefficients
import numpy as np
import draw_eos_uniform as pyeos
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
M_sun_si = const.M_sun.si.value
SLY_CUTOFF_P = 3.9e31 # In SI
sly_polytrope_model = pyeos.eos_polytrope(34.384, 3.005, 2.988, 2.851)
SLY_CUTOFF_EPS = sly_polytrope_model.eval_energy_density(SLY_CUTOFF_P)
SLY_CUTOFF_RHO = sly_polytrope_model.eval_baryon_density(SLY_CUTOFF_P)

# PARAMETRIC EOS
# Model of equation of state prior to sample from: 
# Need to sample gamma0, gamma1, gamma2, gamma3
# From the Carney, Wade, Irwin paper (I don't actually )
gamma0_range = (0.2, 2.0)
gamma1_range = (-1.6, 1.7)
gamma2_range = (-.2, .6)
gamma3_range = (-.02, .02)
parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws")
parser.add_argument("--dir-index", type=int, dest="dir_index")
# need
class eos_spectral:
    def __init__(self,gamma0, gamma1, gamma2, gamma3):
        self.gamma0 = gamma0
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3

        self.eos = [
            self.gamma0,
            self.gamma1, 
            self.gamma2, 
            self.gamma3
        print(gamma0, gamma1, gamma2, gamma3)
        self.logGamma = np.polynomial.polynomial(self.eos)
        self.family = None
        
    # Get the eos family from the paramaters. 
    def eval_P_of_rho(self, rho):
        return rho**(np.exp(self.logGamma))
    def get_eos(self):
        return self.eos
    def solve_tov(self):
        return None
    # Evaluate the eos in terms of epsilon(p)
    def eval_energy_density(self, p):
        integrate()
        
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
    def eval_Gamma(self, p):
        return np.exp(self.gamma0 + self.gamma1 * p + self.gamma2 * p**2 + self.gamma3 * p**3 )
    # Return true if the local speed of sound is larger than the speed of light at the highest pressure allowed for a 
    # certain EOS
    def is_causal(self):
        return True   # Conversion from cgs to SI (leave some leeway like in 1805.11217)
    def is_M_big_enough(self):
        return True

    def is_confined(self, ps):
        return True

             

# We also need gamma1 > gamma2 > gamma3 ? and thermodynamic stability? , so I guess we sample gamma1 first 
# and then constrain the others based on this. This is the (somewhat) uniform prior on the gamma's, I think
# I still need to glue together the 
def get_eos_realization_uniform_spec (gamma0_range = gamma0_range,
                                      gamma1_range= gamma1_range,
                                      gamma2_range=gamma2_range, 
                                      gamma3_range = gamma3_range):
    # There's some problem with configurations not working if the parameters are too close together,
    # so I tried to force them apart without losing too much of the prior
    gamma0 = np.random.uniform(*gamma0_range)
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(*gamma2_range)
    gamma3 = np.random.uniform(*gamma3_range) 
    try:
        print([gamma0, gamma1, gamma2, gamma3])
        return eos_spectral(gamma0, gamma1, gamma2, gamma3)    

    except : 
        # try again :(
        return get_eos_realization_uniform_spec (gamma0_range = gamma0_range,
                                      gamma1_range= gamma1_range,
                                      gamma2_range=gamma2_range, 
                                      gamma3_range = gamma3_range)
        
 # The first prior was not great, most of the EOS's couldn't even 
 # support a 1.4 solar mass neutron star (it really helps to enforce)
 # the criteria ahead of time

# In general I don't think this is the way to go
# def get_eos_realization_improved_para (gamma0_range = gamma0_range,
#                                        gamma1_range= gamma1_range,
#                                        gamma2_range=gamma2_range, 
#                                        gamma3_range = gamma3_range):
#     # Sample around where I know the prior is reasonable
#     eps = .1
#     Cov = np.matrix([[.42,0,0,0],[0,.24,0,0],[0,0,.15,0],[0,0,0,.11]])
#     means = np.array([34.084, 3.205, 2.988, 2.551])
#     samples = np.random.multivariate_normal(means, Cov)
#     # Check if in bounds
#     logp1 = samples[0]
#     gamma3 = samples[1]
#     gamma2 = samples[2]
#     gamma1 = samples[3]
#     g1cond = gamma1_range[0] + eps < gamma1 < gamma2
#     g2cond =  (gamma2_range[0]+eps < gamma2 <gamma3)
#     g3cond =   gamma3 < gamma3_range[1]
#     lpcond = logp1_range[0] < logp1 < logp1_range[1]
#     # Fallback if the criteria aren't satisfied
#     if not (g1cond and g2cond and g3cond and lpcond):
#         print("falling back")
#         return get_eos_realization_uniform_poly()
#     return eos_parametric(logp1, gamma1, gamma2, gamma3)

# Enforce conditions ahead of time
def get_eos_realization_uniform_constrained_spec (gamma0_range = gamma0_range,
                                                  gamma1_range= gamma1_range,
                                                  gamma2_range=gamma2_range,
                                                  gamma3_range = gamma3_range):
    print("gamma0_range is", gamma0_range)
    gamma0 = np.random.uniform(*gamma0_range)
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(gamma2_range[0], gamma2_range[1])
    gamma3 = np.random.uniform(gamma3_range[0], gamma3_range[1])
    
    try :
        this_polytrope = eos_spectral(gamma0, gamma1, gamma2, gamma3)
    except :
        # Try again
        return get_eos_realization_uniform_constrained_spec(gamma0_range = gamma0_range,
                                                            gamma1_range= gamma1_range,
                                                            gamma2_range=gamma2_range,
                                                            gamma3_range = gamma3_range)
    print("Got the polytrope")
    if this_polytrope.is_causal() and this_polytrope.is_M_big_enough():
        return this_polytrope
    else:
        return get_eos_realization_uniform_constrained_spec(gamma0_range = gamma0_range,
                                                            gamma1_range= gamma1_range,
                                      gamma2_range=gamma2_range, gamma3_range = gamma3_range)
    
# Because we have an analytic equation of state, we can compute the derivative dmu/dp 
# analytically.  Therefore we can compute phi analytically (Doesn't seem to actually be necessary)


# Stitch EoS onto the known EoS below nuclear saturation density. 
# Use Sly log(p1) = 34.384, gamma1 = 3.005, gamma2 = 2.988, gamma3 = 2.851
# There's some subtlety here related to where the pressure is known.  Here
# it is given at the dividing line between rho_1 and rho_2,  but we want it
# at the stitching point (I think this is fine though, because we can just 
# evaluate the modelt at the stiching point)

#SLy model
sly_polytrope_model = pyeos.eos_polytrope(34.384, 3.005, 2.988, 2.851)



 
def create_eos_draw_file(name):
    eos_poly = get_eos_realization_uniform_spec(gamma0_range, gamma1_range, gamma2_range, gamma3_range)
    if eos_poly.is_causal() and eos_poly.is_M_big_enough():
    # FIXME WORRY ABOUT CGS VS SI!!!!! (Everything is in SI till the last step :/ ) 
        p_sly = np.geomspace(1.0e30, 3.9e31, 100)
        p_main = np.geomspace (3.9e32, 9.0e35, 800)
        eps_sly = sly_polytrope_model.eval_energy_density(p_sly)
        eps_main = eos_poly.eval_energy_density(p_main)
        rho_b_sly = eos_poly.eval_baryon_density(p_sly)
        rho_b_main = eos_poly.eval_baryon_density(p_main)
        p = np.concatenate([p_sly, p_main])
        eps = np.concatenate([eps_sly, eps_main])
        rho_b = np.concatenate([rho_b_sly, rho_b_main])
        data = np.transpose(np.stack([p/c**2*10 , eps/c**2*10, rho_b/10**3])) # *10 because Everything above is done in SI
        np.savetxt(name,data, header = 'pressurec2, energy_densityc2, baryon_density',
                   fmt='%.10e', delimiter=",")
    else :
        create_eos_draw_file(name)

if __name__ == "__main__":
    args = parser.parse_args()
    num_draws = args.num_draws
    dir_index = args.dir_index
    for i in range(num_draws):
        name = "eos-draw-" + "%06d" %(dir_index*num_draws + i) + ".csv"
        create_eos_draw_file(name)
