import numpy as np
import scipy.interpolate as interp
import lalsimulation as lalsim
import lal
import lalinference as lalinf
import argparse
import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt
import astropy.constants as const
import astropy.units as u
from  scipy.integrate import solve_ivp
import universality.gaussianprocess.utils as gp
import scipy

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
parser.add_argument("--num-draws", type=int, dest="num_draws", default=1)
parser.add_argument("--dir-index", type=int, dest="dir_index", default=0)


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



# START OF SPECTRAL THINGS
# Legacy support, should be removed at some point
p_range = (1e31, 1e37)
p_0 = 3.9e32

# This is the prior used in the introductory paper above
gamma0_range = (0.2, 2.0)
gamma1_range = (-.6, 1.7)
gamma2_range = (-.6, .6)
gamma3_range = (-.02, .02)


# This is the preimage of the prior
# used in https://arxiv.org/pdf/2001.01747.pdf,
# it gets mapped to a relevant range of gammas
r0_range = (-4.37722, 4.91227)
r1_range = (-1.82240, 2.06387)
r2_range = (-.32445, .36469)
r3_range = (-.09529, .11046)


# This maps r's drawn from distributions (such
# as a distributioon with the values above as uniform
# bounds) and returns values of gammma which are more
# likely physical.  This mapping may also be inducing
# correlations in the priors; need to investigate further.
def map_rs_to_gammas(r0, r1, r2, r3):
    S = np.matrix([[.43801, -0.53573, +0.52661, -0.49379],
                   [-0.76705, +0.17169, +0.31255, -0.53336],
                   [+0.45143, 0.67967, -0.19454, -0.54443],
                   [+0.12646, 0.47070, 0.76626, 0.41868]])
    mu_r = np.matrix([[0.89421],[0.33878],[-0.07894],[+0.00393]])
    sigma_r = np.matrix([[0.35700,0,0,0],[0,0.25769,0,0],[0,0,0.05452,0],[0,0,0,0.00312]])
    rs = np.matrix([[r0],[r1], [r2], [r3]])
    return sigma_r * S**(-1) * rs  + mu_r







# This class is meant to hose all of the functions needed to interact with a
# paramaterized eos, without actually exposing the client to any of the lalsimulation
# functions (which can be somewhat volatile and don't come with object orientation)
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
        
############################################################
# Implemented Priors         
############################################################             

# Draws froam a uniform distribution, with the bounds as provided
# in the intro paper, basically useless because of a feature (bug)?
# in lalsim that causes a segfault if parameters produce an EOS
# with insufficient points to do the TOV integration.  So any code that
# runs this function is liable to fail catastrophically
def get_eos_realization_uniform_spec (gamma0_range = gamma0_range,
                                      gamma1_range= gamma1_range,
                                      gamma2_range=gamma2_range, 
                                      gamma3_range = gamma3_range):
    # There's some problem with configurations not working if the parameters are too close together
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
        

# This is something of a wrapper of the lalinference
# function which checks to see if a particular combo
# of parameters will (1) create an EOS which can be used
# to successfully integrate TOV, and (2) create a physically
# plausible solution (perhaps the checks could be relaxed
# to expand the prior but its not really possible to dispense
# with it completely)
def criteria_spec(gamma0, gamma1, gamma2, gamma3):
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
# satisy the criteria, this is also fairly useless because
# most draws are not physical, so it takes an incredibly long
# time to get a sizeable sample of reasonable EOSs
def get_eos_realization_uniform_constrained_spec (gamma0_range = gamma0_range,
                                                  gamma1_range= gamma1_range,
                                                  gamma2_range=gamma2_range,
                                                  gamma3_range = gamma3_range):
    gamma0 = np.random.uniform(*gamma0_range)
    gamma1 = np.random.uniform(*gamma1_range)
    gamma2 = np.random.uniform(gamma2_range[0], gamma2_range[1])
    gamma3 = np.random.uniform(gamma3_range[0], gamma3_range[1])

    try :
        if not criteria_spec(gamma0, gamma1, gamma2, gamma3):
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
    return this_polytrope
    

# Inspired by  https://arxiv.org/pdf/2001.01747.pdf appendix B, see there for help
# This prior uses a separate domain which gets mapped into gamma-space, in doing this
# it targets the most viable regions of parameter space, it is the de facto best choice
# for sampling
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
    failure = False
    try :
        if not criteria_spec(gamma0, gamma1, gamma2, gamma3):
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
    print(gamma0, gamma1, gamma2, gamma3)
    return this_polytrope



# Inspired by  https://arxiv.org/pdf/2001.01747.pdf appendix B, see there for help
def get_eos_realization_mapped_gaussian_constrained_spec (r0_range = r0_range,
                                                 r1_range= r1_range,
                                                 r2_range= r2_range,
                                                 r3_range = r3_range):
    ################################################################
    ranges= [r0_range, r1_range, r2_range, r3_range]
    means = [np.mean(this_range) for this_range in ranges]
    cov = 1/6*np.diag([np.std(this_range) for this_range in ranges])
    [r0, r1, r2, r3] = np.random.multivariate_normal(means, cov) 
    
    

    
    ################################################################
    gammas = map_rs_to_gammas(r0, r1, r2, r3)
    gamma0 = gammas[0,0]
    gamma1 = gammas[1,0]
    gamma2 = gammas[2,0]
    gamma3 = gammas[3,0]
    failure = False
    try :
        if not criteria_spec(gamma0, gamma1, gamma2, gamma3):
            failure = True
        this_polytrope = eos_spectral(gamma0, gamma1, gamma2, gamma3)
        if 3 < this_polytrope.get_max_M() or 1.9 > this_polytrope.get_max_M():
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
    # This is a lot of printing, but makes it possible to diagnose the prior more easily
    print(gamma0, gamma1, gamma2, gamma3)
    return this_polytrope









# Stitch EoS onto the known EoS below nuclear saturation density. 
# Use Sly log(p1) = 34.384, gamma1 = 3.005, gamma2 = 2.988, gamma3 = 2.851
# There's some subtlety here related to where the pressure is known.  Here
# it is given at the dividing line between rho_1 and rho_2,  but we want it
# at the stitching point (I think this is fine though, because we can just 
# evaluate the modelt at the stiching point)



sly_polytrope_model = eos_polytrope(34.384, 3.005, 2.988, 2.851)

def get_table_small_vals(file_name, rho_thresh):
    data = np.loadtxt(file_name, delimiter=",", skiprows=1)
    threshold_index = np.min(np.where(data[:,2] > rho_thresh))
    return data[:threshold_index, :]

def sample_table(eos_dir, num_dirs, draws_per_dir, rho_thresh = 4.8e14):
    fast_vary = np.random.randint(draws_per_dir)
    slow_vary = np.random.randint(num_dirs)
    dir_name = "DRAWmod" + str(draws_per_dir)+"-" + str(slow_vary).zfill(6)
    file_name = "/eos-draw-"+ str(draws_per_dir*slow_vary + fast_vary).zfill(6) + ".csv"
    file_path = eos_dir + dir_name + file_name
    data = get_table_small_vals(file_path, rho_thresh)
    # Actually switch to SI to avoid spaghetti code
    p_small = data[:, 0]/10*c**2
    eps_small = data[:,1]/10*c**2
    rho_b_small = data[:,2]*10**3
    return p_small, eps_small, rho_b_small
def this_cs2(eps):
    return .75
# Get the extremal baryon density and pressure from known speed of sound
def solve_for_quants(eps_vals, p0, cs_2): 
    def dZ(eps, Z):
        deriv = np.zeros((1))
        p = Z[0]
  
        # p deriv, dp
        deriv[0] = cs_2(eps)

        return deriv
    eps0 = min(eps_vals)-10**9
    epsf = max(eps_vals)+10**9
    Z0 = np.array([p0])
    sol= solve_ivp(dZ, (eps0,epsf), Z0, method="Radau", dense_output="True", t_eval=eps_vals)
    eps_main = sol.t
    p_main = sol.y[0,:]
    
    print(eps_main.size)
    print(p_main.size)


    return eps_main, p_main

def e_p2rho(epsc2, pc2, rhob0):
    dlogrhodeps = 1/(epsc2 + pc2)
    DeltaLogrho = scipy.integrate.cumulative_trapezoid(dlogrhodeps, epsc2, initial=0)
    print("size of output of cumulative_trapezoid is", DeltaLogrho.shape)
    return rhob0*np.exp(DeltaLogrho)
# Just a wrapper to keep things more general overall sly has nothing to do with this
def sample_small_data_nonparametric(sly=False):
    return sample_table("/home/philippe.landry/nseos/eos/gp/mrgagn/", 2000, 1000, rho_thresh=4.8e14)
def sample_small_data_piecewise(sly=False):
    eos_poly = None
    if sly :
        eos_poly = sly_polytrope_model
    
        # Want denser spacing at high values
        p_small = np.geomspace(2e33, 1.0e12, 750)
        # Orient in the correct direction
        p_small = np.flip(p_small)
        eps_small=  eos_poly.eval_energy_density(p_small)
        rho_b_small = eos_poly.eval_baryon_density(p_small)
        return p_small, eps_small, rho_b_small
    else:
        return sample_table("/home/isaac.legred/parametric-eos-priors/production_eos_draw_piecewise/", 2000, 100)
def sample_small_data_spectral(sly=False):
    eos_spec = None
    # This is not good
    if sly :
        eos_spec = sly_polytrope_model
        # Want denser spacing at high values
        p_small = np.geomspace(2e33, 1.0e12, 750)
        # Orient in the correct direction
        p_small = np.flip(p_small)
        # Walks like a duck, talks like a duck
        eps_small=  eos_spec.eval_energy_density(p_small)
        rho_b_small = eos_spec.eval_baryon_density(p_small)
        return p_small, eps_small, rho_b_small
    else:
        return sample_table("/home/isaac.legred/parametric-eos-priors/production_eos_draw_spectral/", 1000, 200)
    
def create_eos_draw_file(name, small_data_getter = sample_small_data_piecewise,  sly=False):
    p_small, eps_small, rho_b_small = small_data_getter(sly)
    deltaeps=0
    eps_transp = (1+ deltaeps) *eps_small[-1]
    # get the jump in rho across the transition
    logrho_trans = (np.log(eps_transp + p_small[-1]) - np.log(eps_small[-1] + p_small[-1]))
    eps_main = np.geomspace(eps_transp, 10.0**37, 750)
    print("eps_main[0] = ", eps_main[0])
    eps_main, p_main = solve_for_quants(eps_main, p_small[-1], this_cs2)
    #p_of_rhob = scipy.interpolate.interp1d(rho_b_small, p_small, kind="cubic")
   
    p_mainc2 = p_main/c**2*10
    eps_mainc2 = eps_main/c**2*10
    pc2 = np.concatenate([p_small, p_main])*10/c**2 # To cgs
    epsc2 = np.concatenate([eps_small, eps_main])*10/c**2 # To cgs
   


    rho_b_main=e_p2rho(eps_mainc2,  p_mainc2, rho_b_small[-1]/1000 *(np.exp(logrho_trans))) # in cgs
    rho_b = np.concatenate([rho_b_small/10**3, rho_b_main]) # in cgs
    print(rho_b.shape)
    data = np.transpose(np.stack([pc2 , epsc2, rho_b]))
    np.savetxt(name,data, header = 'pressurec2,energy_densityc2,baryon_density',
               fmt='%.10e', delimiter=",", comments="")

def create_sly_draw_file(name):
    create_eos_draw_file(name, sly=True)



if __name__ == "__main__":
    args = parser.parse_args()
    num_draws = args.num_draws
    dir_index = args.dir_index
    for i in range(num_draws):
        name = "eos-draw-" + "%06d" %(dir_index*num_draws + i) + ".csv"
        if i!=0:
            create_eos_draw_file(name, small_data_getter=sample_small_data_nonparametric)
        else:
            create_sly_draw_file(name)
