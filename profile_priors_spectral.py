# I guess this should probably be refactored so that all of the writing functions are independent of the sampling 
# functions
import draw_eos_spectral  as py_spec_eos
import draw_eos_uniform as pyeos
import lalsimulation as lalsim
import lal
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
# Number of samples to draw
N = 20
import astropy.constants as const
Gsi = const.G.si.value
csi  = const.c.si.value
ccgs = const.c.cgs.value
eos_list = []
for i in range(N):
    eos = py_spec_eos.get_eos_realization_uniform_constrained_spec()  #Use default parameter range
    eos.get_fam()
    eos_list.append(eos) 
    print("eos is causal : ", eos.is_causal())
    print("eos is confined : ", eos.is_confined(np.linspace(*py_spec_eos.p_range, 100)))
    print("i is ", i)
    
max_masses = np.array([lalsim.SimNeutronStarMaximumMass(eos.get_fam())/lal.MSUN_SI for eos in eos_list if eos.get_fam() is not None])


M1p4s = lal.MSUN_SI*1.4
R1p4s =[] 
nuc_rho_cgs = 2.8e14 # This is in cgs
nuc_rho_si = 2.8e17 # this is in SI

# I can't believe there's no function to get pressure as a function of density.
def find_pressure_of_density(rho, eos_wrapper):
    def rootable(logp):
        p = np.exp(logp)
        return  np.arcsinh(eos_wrapper.eval_baryon_density(p) - rho)
    return np.exp((scipy.optimize.root_scalar(rootable, bracket=(np.log(2.0e27), np.log(1.0e37))).root))
P2nucs = []

# Testing ###############
p = 1e37
eos_0 = eos_list[0]
print("p2nuc missed by", eos_0.eval_baryon_density(p) - 2*nuc_rho_si)

#########################
for  n, eos in enumerate (eos_list):
    # There has to be a better way to do this in the long run
    #P2nucs.append(np.log10((find_pressure_of_density(nuc_rho_si, eos))/ccgs**2*10)) # SI -> cgs is x 10
    if eos is not None:
        if max_masses[n] > 1.43:
            try:
                R1p4s.append(lalsim.SimNeutronStarRadius(M1p4s, eos.family)/1000.)
                P2nucs.append(np.log10((find_pressure_of_density(2 * nuc_rho_si,
                                                                 eos))/ccgs**2*10))
                                                          # SI -> cgs is x 10 In cgs)
            except:
                pass
        else:
            print("failed to produce a reasonable max mass")
#Loves = np.array([lalsim.SimNeutronStarLoveNumberK2(M1p4s, eos.family) for eos in eos_list])

 #Mass plot
plt.hist(max_masses, bins=30)
plt.xlabel("max_mass")
plt.ylabel("frequency")
plt.savefig("better_prior_masses.png")

# R_1.4
plt.figure()
plt.hist(R1p4s, bins=30)
plt.xlabel(r"$R_{1.4}$")
plt.ylabel("frequency")
plt.savefig("better_prior_radii.png")

print(P2nucs)
# Pressure(2 rho_nuc)
plt.figure()
plt.hist(P2nucs, bins=30)
plt.xlabel(r"$\log_{10}(P(2\rho_{nuc})/c^2$)")
plt.ylabel("frequency")
plt.savefig("better_prior_P2nuc.png")
plt.show()

