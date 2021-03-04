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
# Need to sample gamma0, gamma1, gamma2, gamma3
# From the Carney, Wade, Irwin paper (I don't actually )

p_range = (1e31, 1e37)
p_0 = 3.9e32

gamma0_range = (0.2, 2.0)
gamma1_range = (-.6, 1.7)
gamma2_range = (-.6, .6)
gamma3_range = (-.02, .02)

r0_range = (-4.37722, 4.91227)
r1_range = (-1.82240, 2.06387)
r2_range = (-.32445, .36469)
r3_range = (-.09529, .11046)

parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws")
parser.add_argument("--dir-index", type=int, dest="dir_index")
# need

def map_rs_to_gammas(r0, r1, r2, r3):
    S = np.matrix([[.43801, -0.53573, +0.52661, -0.49379],
                   [-0.76705, +0.17169, +0.31255, -0.53336],
                   [+0.45143, 0.67967, -0.19454, -0.54443],
                   [+0.12646, 0.47070, 0.76626, 0.41868]])
    mu_r = np.matrix([[0.89421],[0.33878],[-0.07894],[+0.00393]])
    sigma_r = np.matrix([[0.35700,0,0,0],[0,0.25769,0,0],[0,0,0.05452,0],[0,0,0,0.00312]])
    rs = np.matrix([[r0],[r1], [r2], [r3]])
    return sigma_r * S**(-1) * rs  + mu_r
class eos_spectral:
    def __init__(self,gamma0, gamma1, gamma2, gamma3):
        self.gamma0 = gamma0
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3

        self.eos = lalsim.SimNeutronStarEOS4ParameterSpectralDecomposition(
            self.gamma0,
            self.gamma1, 
            self.gamma2, 
            self.gamma3)
        print(gamma0, gamma1, gamma2, gamma3)
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

    def is_confined(self, ps):
        if (.6 < self.eval_Gamma(ps).all() < 4.5):
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
# Enforce conditions ahead of time
def criteria(gamma0, gamma1, gamma2, gamma3):
    vars = lalinf.Variables()
    no_vary = lalinf.lalinference.LALINFERENCE_PARAM_FIXED
    lalinf.lalinference.AddREAL8Variable(vars, "SDgamma0", gamma0, no_vary )
    lalinf.lalinference.AddREAL8Variable(vars, "SDgamma1", gamma1, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "SDgamma2",  gamma2, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "SDgamma3",  gamma3, no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "mass1",  1.4 , no_vary)
    lalinf.lalinference.AddREAL8Variable(vars, "mass2",  1.4 , no_vary)
    
    a = lal.CreateStringVector("Hi")
    process_ptable = lalinf.ParseCommandLineStringVector(a)
    success_param = lalinf.EOSPhysicalCheck(vars, process_ptable)
    if success_param == 0:
        return True
    else :
        return False
# Get an EOS sampled from a uniform prior that is supposed to 
# satisy the criteria
def get_eos_realization_uniform_constrained_spec (gamma0_range = gamma0_range,
                                                  gamma1_range= gamma1_range,
                                                  gamma2_range=gamma2_range,
                                                  gamma3_range = gamma3_range):
    gamma0 = np.random.uniform(*gamma0_range)
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(gamma2_range[0], gamma2_range[1])
    gamma3 = np.random.uniform(gamma3_range[0], gamma3_range[1])

    try :
        if not criteria(gamma0, gamma1, gamma2, gamma3):
            raise ValueError("Did not meet the required Gamma Criteria") 
        this_polytrope = eos_spectral(gamma0, gamma1, gamma2, gamma3)
        # In principle the maximum mass below 1.9 should be fine,
        # but there is no way to know 
        if 3 < this_polytrope.get_max_M() or 1.7 > this_polytrope.get_max_M():
            print("throwing an exception")
            raise ValueError("M is not in the right range")
        print(":)")
    except :
        print(":)")
        # Try again
        return get_eos_realization_uniform_constrained_spec(gamma0_range = gamma0_range,
                                                            gamma1_range= gamma1_range,
                                                            gamma2_range=gamma2_range,
                                                            gamma3_range = gamma3_range)
    
    # Might want this condition at some point
    #and this_polytrope.is_M_big_enough()
    return this_polytrope
    
# Because we have an analytic equation of state, we can compute the derivative dmu/dp 
# analytically.  Therefore we can compute phi analytically (Doesn't seem to actually be necessary)

# Inspired by  https://arxiv.org/pdf/2001.01747.pdf appendix B, see there for help
def get_eos_realization_mapped_constrained_spec (r0_range = r0_range,
                                                 r1_range= r1_range,
                                                 r2_range= r2_range,
                                                 r3_range = r3_range):
    r0 = np.random.uniform(*r0_range)
    r1 = np.random.uniform(*r1_range)
    r2 = np.random.uniform(*r2_range)
    r3 = np.random.uniform(*r3_range)

    gammas = map_rs_to_gammas(r0, r1, r2, r3)
    gamma0 = gammas[0,0]
    gamma1 = gammas[1,0]
    gamma2 = gammas[2,0]
    gamma3 = gammas[3,0]
    print(gamma0, gamma1, gamma2, gamma3)
    failure = False
    try :
        if not criteria(gamma0, gamma1, gamma2, gamma3):
            failure = True
        this_polytrope = eos_spectral(gamma0, gamma1, gamma2, gamma3)
        if 3 < this_polytrope.get_max_M() or 1.7 > this_polytrope.get_max_M():
            failure = True
        
    except :
        # Try again
        return get_eos_realization_mapped_constrained_spec(r0_range = r0_range,
                                                            r1_range= r1_range,
                                                            r2_range=r2_range,
                                                            r3_range = r3_range)
    if failure:
         return get_eos_realization_mapped_constrained_spec(r0_range = r0_range,
                                                            r1_range= r1_range,
                                                            r2_range= r2_range,
                                                            r3_range = r3_range)
    return this_polytrope



# Inspired by  https://arxiv.org/pdf/2001.01747.pdf appendix B, see there for help
def get_eos_realization_mapped_gaussian_constrained_spec (r0_range = r0_range,
                                                 r1_range= r1_range,
                                                 r2_range= r2_range,
                                                 r3_range = r3_range):
    ################################################################
    ranges= [r0_range, r1_range, r2_range, r3_range]
    # Do something cleverer here
    means = [np.mean(this_range) for this_range in ranges]
    cov = 1/9*np.diag([np.std(this_range) for this_range in ranges])
    [r0, r1, r2, r3] = np.random.multivariate_normal(means, cov) 
    
    

    
    ################################################################
    gammas = map_rs_to_gammas(r0, r1, r2, r3)
    gamma0 = gammas[0,0]
    gamma1 = gammas[1,0]
    gamma2 = gammas[2,0]
    gamma3 = gammas[3,0]
    print(gamma0, gamma1, gamma2, gamma3)
    failure = False
    try :
        if not criteria(gamma0, gamma1, gamma2, gamma3):
            failure = True
        this_polytrope = eos_spectral(gamma0, gamma1, gamma2, gamma3)
        if 3 < this_polytrope.get_max_M() or 1.7 > this_polytrope.get_max_M():
            failure = True
    except :
        # Try again
        return get_eos_realization_mapped_gaussian_constrained_spec(r0_range = r0_range,
                                                            r1_range= r1_range,
                                                            r2_range=r2_range,
                                                            r3_range = r3_range)
    if failure:
        return get_eos_realization_mapped_gaussian_constrained_spec(r0_range = r0_range,
                                                            r1_range= r1_range,
                                                            r2_range= r2_range,
                                                            r3_range = r3_range)
    return this_polytrope









# Stitch EoS onto the known EoS below nuclear saturation density. 
# Use Sly log(p1) = 34.384, gamma1 = 3.005, gamma2 = 2.988, gamma3 = 2.851
# There's some subtlety here related to where the pressure is known.  Here
# it is given at the dividing line between rho_1 and rho_2,  but we want it
# at the stitching point (I think this is fine though, because we can just 
# evaluate the modelt at the stiching point)

#SLy model
sly_polytrope_model = pyeos.eos_polytrope(34.384, 3.005, 2.988, 2.851)



 
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
