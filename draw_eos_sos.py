#Test Sampling eos's from spectral coefficients
import numpy as np
import draw_eos_piecewise as pyeos
import scipy.interpolate as interp
import lalsimulation as lalsim
import lalinference as lalinf
import lal
import argparse
from matplotlib import pyplot as plt
import astropy.constants as const
import astropy.units as u


c = const.c.cgs.value
# Characteristic refinement number (probably should be changed on a case by case basis)
N = 100
M_sun_si = const.M_sun.si.value

# PARAMETRIC EOS
# Model of equation of state prior to sample from: 


#https://arxiv.org/pdf/1812.08188.pdf


# Legacy support, should be removed at some point
p_range = (1e31, 1e37)
p_0 = 3.9e32


# What to put here?
rho_0 = 2.4e14 # g/cm**3

a1_range=
a2_range=
a3_range=
a4_range=
a5_range=
#a6 is special
#a6_range=

# This can be called as a script, in that case it produces a single "draw file" which contains
# a tabulated eos of (pressure, energy density, baryon density)
parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws")
parser.add_argument("--dir-index", type=int, dest="dir_index")
def get_cs2(a1, a2, a3, a4, a5, a6):
    fun = lambda x : a1*np.exp(-1/2*(x-a2)**2/a3**2) + a6 +(1/3 - a6)/(1 + e**(-a5(x-a4)))
    return fun
def tabulate_values(eps_min, eps_max, cs_2):
    # Find low density eos
    # they glue to a particular eos but I think it's better
    # to glue to SLy for consistency.  
    # Somehow need to know the value of a before this point?
# This class is meant to hose all of the functions needed to interact with a
# paramaterized eos, without actually exposing the client to any of the lalsimulation
# functions (which can be somewhat volatile and don't come with object orientation)
sly_polytrope_model = pyeos.eos_polytrope(34.384, 3.005, 2.988, 2.851)
sly_matching_eps = py_eos.
class eos_speed_of_sound:
    def __init__(self, a1, a2, a3, a4, a5, a6):

        self.x = #eps/(m_N n_0)
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5
        self.a6 = a6

        self.cs2 = get_cs2(self.a1, self.a2, self.a3, self.a4, self.a5 self.a6)
        p, eps, rho_b = tabulate_values(eps_min, eps_max, cs2)
        
        self.family = lalsim.CreateSimNeutronStarFamily(self.eos)
        
    # Get the eos family from the paramaters. 
    def get_eos(self):
        return self.eos
    def get_fam(self):
        return self.family
    # Evaluate the eos in terms of epsilon(p)
    def eval_energy_density(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):    
            eps = np.zeros(len(p))  
            for i, pres in enumerate(p):
                eps[i] = lalsim.SimNeutronStarEOSEnergyDensityOfPressure(pres,self.eos)    
        else:
            eps = lalsim.SimNeutronStarEOSEnergyDensityOfPressure(p, self.eos)
        return eps
    # Evaluate the phi parameter as used in the non-parametric papers,
    # Not currently used
    def eval_phi(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):    
            eps = np.zeros(len(p))  
            for i, pres in enumerate(p):
                eps[i] = lalsim.SimNeutronStarEOSEnergyDensityDerivOfPressure(pres,self.eos)    
        else:
             eps = lalsim.SimNeutronStarEOSEnergyDensityDerivOfPressure(p, self.eos)
        return eps
    # Evaluate the baryon density at a particular pressure
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
    # Evluate the speed of sound at a particular pressure
    def eval_speed_of_sound(self, p):
        if isinstance(p, list) or isinstance(p, np.ndarray):
            cs = np.zeros(len(p))
            for i, pres in enumerate(p):
                try:
                    h = lalsim.SimNeutronStarEOSPseudoEnthalpyOfPressure(pres, self.eos)
                    cs[i] = lalsim.SimNeutronStarEOSSpeedOfSound(h, self.eos)
                except:
                    print(pres, "failed to produce a valid sound speed")
                    break
        else:
            cs  = lalsim.SimNeutronStarEOSSpeedOfSound(p, self.eos)
        return cs
    #Evaluate the exponent polytope (I don't really know if this works, wouldn't recommend using)
    def eval_Gamma(self, p):
        x = np.log( p/p_0)
        return np.exp(self.gamma0 + self.gamma1 * x + self.gamma2 * x**2 + self.gamma3 * x**3 )
    # Return true if the local speed of sound is larger than the speed of light at the highest pressure allowed for a 
    # certain EOS
    def is_causal(self, ps):
        c_si = c/100 # The speed of sound in SI
        cs = self.eval_speed_of_sound(ps)
        cs_max = max(cs)
        print("cs_max is", cs_max)
        return cs_max < c_si*1.1
    def get_max_M(self):
        return lalsim.SimNeutronStarMaximumMass(self.family)/lal.MSUN_SI 
    # This function claims to check to see if the adiabatic exponent is bounded in the
    # range [.6, 4.5], it only works as well as the function which evaluates the gamma
    def is_confined(self, ps):
        if (.6 < self.eval_Gamma(ps).all() < 4.5):
            return True
        



 
def create_eos_draw_file(name):
    print(name)
    eos_poly = get_eos_realization_mapped_gaussian_constrained_spec(r0_range, r1_range, r2_range, r3_range)
    if True:
    # FIXME WORRY ABOUT CGS VS SI!!!!! (Everything is in SI till the last step :/ ) 
        p_small = np.linspace(1.0e12, 1.3e30, 300)
        p_main = np.linspace (1.3e30, 9.0e36, 500)
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
    else :
        create_eos_draw_file(name)

if __name__ == "__main__":
    args = parser.parse_args()
    num_draws = args.num_draws
    dir_index = args.dir_index
    for i in range(num_draws):
        name = "eos-draw-" + "%06d" %(dir_index*num_draws + i) + ".csv"
        create_eos_draw_file(name)
