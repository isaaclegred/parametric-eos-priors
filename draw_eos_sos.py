#Sampling EoSs from a parametrization of the speed of sound.
import numpy as np
import draw_eos_piecewise as pyeos
import scipy.interpolate as interp
import scipy.optimize as optimize
import scipy.integrate as integrate
import scipy
import lalsimulation as lalsim
import lalinference as lalinf
import lal
import argparse
import matplotlib as mpl
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



# We want to transition to the model slightly below
rho_0 = 2.8e14 # g/cm**3
rho_0_si = 2.8e17 # kg/m**3
rho_small = rho_0/2
c2si = (2.998e8)**2

a1_range=(.1, 1.5)
a2_range=(1.5, 12)
a3_range=(.05, 2)
a4_range=(1.5, 37)
a5_range=(.1, 1)
#a6 is special
#a6_range=

# This can be called as a script, in that case it produces a single "draw file" which contains
# a tabulated eos of (pressure, energy density, baryon density)
parser = argparse.ArgumentParser(description='Get the number of draws needed, could be expanded')
parser.add_argument("--num-draws", type=int, dest="num_draws", default=1)
parser.add_argument("--dir-index", type=int, dest="dir_index", default=0)
parser.add_argument("--prior-tag", type=str, dest="prior_tag", default="uniform")

sly_polytrope_model = pyeos.eos_polytrope(34.384, 3.005, 2.988, 2.851)
sly_p_1 = 2e32 #(SI, i.e. Pa)  This is an arbitrary cutoff point
sly_eps_1 = sly_polytrope_model.eval_energy_density(sly_p_1) #SI J/cm**3
sly_rho_1 = sly_polytrope_model.eval_baryon_density(sly_p_1) #SI g/cm**3

# Eval f1(val) if val > thresh and f2(val) otherwise
# doesn't work with arrays
def function_switch(high_val_func, low_val_func, thresh, val):
    if val > thres:
        return high_val_func(val)
    else:
        return low_val_func(val)

# Return a function which gives the speed of sound at all points
def get_cs2c2(a1, a2, a3, a4, a5, a6):
    print("something broken in get_cs2c2?", a1, a2, a3, a4, a5, a6)
    fun = lambda x : a1*np.exp(-1/2*(x-a2)**2/a3**2) + a6 +(1/3 - a6)/(1 + np.exp(-a5*(x-a4)))
    return fun
def tabulate_values(eps_min, eps_max, cs_2c2, p_min, rho_min):
    # Find low density eos
    # they glue to a particular eos but I think it's better
    # to glue to SLy for consistency.  
    # p_min should be decided upon ahead of time, and (eps_min, p_min)
    # should be a point on the SLy eps-p curve which is low enough
    # to be considered "known"
    eps_vals = np.geomspace(eps_min, eps_max, 500)
    def dp_and_rho(eps, p_and_rho):
        p = p_and_rho[0]
        rho = p_and_rho[1]
        # Define the pressure and density differentials
        dp = cs_2c2(eps/rho_0/c2si)
        drho = rho/(eps + p) # Don't think too hard about the units
        return np.array([dp, drho])
    tabulated_eos = integrate.solve_ivp(dp_and_rho, ((1-10**-10)*eps_min, (1+10**-10)*eps_max), np.array([p_min, rho_min]), t_eval=eps_vals)
    eps = tabulated_eos.t
    p = tabulated_eos.y[0,:]
    rho =tabulated_eos.y[1,:]
    # Note the order of the returns
    return p, eps, rho

# An EoS model based on the speed of sound parametrization

class eos_speed_of_sound:
    def compute_a6(self):
        # In the paper this is made to match some EFT, but here we want to match it to SLy to make
        # it consistent at low denisites with the other EoSs
        # Shoot for it?
        to_match  = self.sly_model.eval_speed_of_sound(sly_p_1)**2/c2si
        eps_match = self.sly_model.eval_energy_density(sly_p_1)
        def diff(self, a_6_guess):
            cs2_guess = self.construct_cs2_helper(a_6_guess)
            print("to_match is", to_match)
            print("eps_match is", eps_match/(rho_0_si*c2si))
            print("guess is",cs2_guess(eps_match/(rho_0_si*c2si)))
            error = to_match - cs2_guess(eps_match/(c2si*rho_0_si))
            print("error is", error)
            return error
        result = optimize.root_scalar(lambda a_6_guess : diff(self, a_6_guess), x0=0, x1=1, fprime = lambda x : -1)
        print("root is", result.root)
        return result.root

    def __init__(self, a1, a2, a3, a4, a5):

        self.x =np.linspace(0,16,1000) #eps/(m_N n_0)
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5
        self.sly_model= sly_polytrope_model
        self.eps_of_p = None
        self.rho_of_p = None
        # The lower point is  arbitrary
        self.p_small = np.geomspace(10**6, sly_p_1, 200)
        self.eps_small = sly_polytrope_model.eval_energy_density(self.p_small)
        self.rho_small = sly_polytrope_model.eval_baryon_density(self.p_small)
        self.p_main = None
        self.eps_main = None
        self.rho_main = None
        # need to find a6 by gluing, should have a procedure to do this right
        # off the bat i think (edit : implemented, need to check if works)
        self.a6  = self.compute_a6()
        self.cs2 = self.construct_cs2_helper(self.a6)
        eps_min = sly_eps_1
        eps_max = 7e20*c2si # in SI, kg/m^3, relatively arbitrary
        self.p_main, self.eps_main, self.rho_main = tabulate_values(eps_min, eps_max, self.cs2, sly_p_1, sly_rho_1)
        

    # Evaluate the eos in terms of epsilon(p)
    def eval_energy_density(self, p, use_low_without_eval=False):
        # if you know you're using the low-p definition ahead of time
        if use_low_without_eval:
            return sly_polytrope_model.eval_energy_density(p1)
        if self.eps_of_p == None:
            #Interpolate a function eps(p)
            self.eps_of_p = scipy.interpolate.interp1d(np.concatenate([self.p_small, self.p_main]), np.concatenate([self.eps_small, self.eps_main]), kind='cubic')
        return self.eps_of_p(p)
    # Evaluate the baryon density at a particular pressure
    def eval_baryon_density(self, p):
        #interpolate a function rho(p)
        if self.rho_of_p is None:
            self.rho_of_p = scipy.interpolate.interp1d(np.concatenate([self.p_small, self.p_main]), np.concatenate([self.rho_small, self.rho_main]), kind='cubic')
        return self.rho_of_p(p)
    # Evluate the speed of sound at a particular pressure
    def eval_speed_of_sound(self, p):
        return self.cs2(self.eval_energy_density(p)/rho_0)
    # thin wrapper around function that actually computes the thing (rename?)
    def construct_cs2_helper(self, a_6):
        a1 = self.a1
        a2 = self.a2
        a3 = self.a3
        a4 = self.a4
        a5 = self.a5
        a6 = a_6
        print("a_6 in construct_cs2 is", a_6)
        return get_cs2c2(a1, a2, a3, a4, a5, a6)
            



def get_eos_realization_sos(a1_range=a1_range, a2_range=a2_range, a3_range=a3_range, a4_range=a4_range, a5_range=a5_range):
    a1=np.random.uniform(*a1_range)
    a2=np.random.uniform(*a2_range)
    a3=np.random.uniform(*a3_range)
    a4=np.random.uniform(*a4_range)
    a5=np.random.uniform(*a5_range)
    return eos_speed_of_sound(a1, a2, a3, a4, a5)
def create_eos_draw_file(name):
    print(name)
    eos_poly = get_eos_realization_sos(a1_range, a2_range, a3_range, a4_range, a5_range)
    if True:
    # FIXME WORRY ABOUT CGS VS SI!!!!! (Everything is in SI till the last step  ) 
        p_small = np.geomspace(1.2e10, 1.3e30, 300)
        p_main = np.geomspace (1.4e30, 9.0e36, 500)
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
