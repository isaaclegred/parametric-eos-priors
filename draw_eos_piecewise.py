# Test Sampling eos's from parametric values
import numpy as np
import scipy.interpolate as interp
import lalsimulation as lalsim
import lal
import lalinference as lalinf
import argparse
from matplotlib import pyplot as plt
import astropy.constants as const
import astropy.units as u


c = const.c.cgs.value
# Characteristic refinement number (probably should be changed on a case by case basis)
N = 200
M_sun_si = const.M_sun.si.value

# Model of equation of state prior to sample from: 
# Need to sample gamma1, gamma2, gamma3, and logp1
# From the Lackey and Wade paper (I don't actually )
logp1_range = (33.6, 35.4)
gamma1_range  = (1.9, 4.5)
gamma2_range = (1.1, 4.5)
gamma3_range = (1.1, 4.5)

parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws")
parser.add_argument("--dir-index", type=int, dest="dir_index")
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
    # Return true if the local speed of sound is larger than the speed of light at the highest pressure allowed for a 
    # certain EOS
    def is_causal(self):
        p_max = lalsim.SimNeutronStarEOSMaxPressure(self.eos)
        c_s_max= lalsim.SimNeutronStarEOSSpeedOfSound(
                    lalsim.SimNeutronStarEOSPseudoEnthalpyOfPressure(p_max,self.eos), self.eos)
        return c_s_max < c/100*1.1   # Conversion from cgs to SI (leave some leeway like in 1805.11217)
    def is_M_big_enough(self):
        m_max = lalsim.SimNeutronStarMaximumMass(self.family)
        return m_max > 1.76 * M_sun_si
        

             

# We also need gamma1 > gamma2 > gamma3 ? and thermodynamic stability? , so I guess we sample gamma1 first 
# and then constrain the others based on this. This is the (somewhat) uniform prior on the gamma's, I think
# I still need to glue together the 
def get_eos_realization_uniform_poly (logp1_range = logp1_range, gamma1_range= gamma1_range, gamma2_range=gamma2_range, 
                                        gamma3_range = gamma3_range):
    # There's some problem with configurations not working if the parameters are too close together,
    # so I tried to force them apart without losing too much of the prior
    eps = .1
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(gamma2_range[0]+eps, gamma1 - eps)
    gamma3 = np.random.uniform(gamma3_range[0], gamma2 - eps) 
    logp1 = np.random.uniform(*logp1_range)
    return eos_polytrope(logp1, gamma1, gamma2, gamma3)    
 

# Enforce conditions ahead of time
def criteria(logp1, gamma1, gamma2, gamma3):
    vars = lalinf.Variables()
    no_vary = lalinf.lalinference.LALINFERENCE_PARAM_FIXED
    lalinf.lalinference.AddREAL8Variable(vars, "logp1", logp1, no_vary )
    lalinf.lalinference.AddREAL8Variable(vars, "gamma1", gamma1, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "gamma2",  gamma2, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "gamma3",  gamma3, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "mass1",  1.4 , no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "mass2",  1.4 , no_vary)
    
    a = lal.CreateStringVector("Hi")
    process_ptable = lalinf.ParseCommandLineStringVector(a)
    success_param = lalinf.EOSPhysicalCheck(vars, process_ptable)
    if success_param == 0:
        return True
    else :
        return False


def get_eos_realization_uniform_constrained_poly (logp1_range = logp1_range,
                                                  gamma1_range= gamma1_range,
                                                  gamma2_range=gamma2_range,
                                                  gamma3_range = gamma3_range):
    # There's some problem with configurations not working if the parameters are too close together,
    # so I tried to force them apart without losing too much of the prior
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(*gamma2_range)
    gamma3 = np.random.uniform(*gamma3_range)
    logp1 = np.random.uniform(*logp1_range)
    this_polytrope = eos_polytrope(logp1, gamma1, gamma2, gamma3)
    if criteria(logp1, gamma1, gamma2, gamma3):
        return this_polytrope
    else:
        return get_eos_realization_uniform_constrained_poly(logp1_range = logp1_range,
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
sly_polytrope_model = eos_polytrope(34.384, 3.005, 2.988, 2.851)



 
def create_eos_draw_file(name):
    eos_poly = get_eos_realization_uniform_constrained_poly(logp1_range, gamma1_range, gamma2_range, gamma3_range)
    # FIXME WORRY ABOUT CGS VS SI!!!!! (UPDATE : Everything is in SI till the last step :/ ) 
    p_small = np.linspace(1.0e12, 1.3e30, 600)
    p_main = np.geomspace (1.3e30, 9.0e36, 700)
    eps_small=  eos_poly.eval_energy_density(p_small)
    eps_main = eos_poly.eval_energy_density(p_main)
    rho_b_small = eos_poly.eval_baryon_density(p_small)
    rho_b_main = eos_poly.eval_baryon_density(p_main)
    p = np.concatenate([p_small, p_main])
    eps = np.concatenate([eps_small, eps_main])
    rho_b = np.concatenate([rho_b_small, rho_b_main])
    data = np.transpose(np.stack([p/c**2*10 , eps/c**2*10, rho_b/10**3])) # *10 because Everything above is done in SI
    np.savetxt(name,data, header = 'pressurec2,energy_densityc2,baryon_density',
               fmt='%.10e', delimiter=",", comments="")
    

if __name__ == "__main__":
    args = parser.parse_args()
    num_draws = args.num_draws
    dir_index = args.dir_index
    for i in range(num_draws):
        name = "eos-draw-" + "%06d" %(dir_index*num_draws + i) + ".csv"
        create_eos_draw_file(name)
