# I guess this should probably be refactored so that all of the writing functions are independent of the sampling 
# functions
import draw_eos_uniform as pyeos
import lalsimulation as lalsim
import lal
import numpy as np
import matplotlib.pyplot as plt
# Number of samples to draw
N = 500
import astropy.constants as const
Gsi = const.G.si.value
csi  = const.c.si.value
ccgs = const.c.cgs.value
eos_list = []
for i in range(N):
    eos_list.append(pyeos.get_eos_realization_improved_poly()) #Use default parameter range
max_masses = np.array([lalsim.SimNeutronStarMaximumMass(eos.family)/lal.MSUN_SI for eos in eos_list])
M1p4s = lal.MSUN_SI*1.4
R1p4s =[] 
nuc_rho_cgs = 2.8e14 # This is in cgs
nuc_rho_si = 2.8e17*csi**2 # this is in SI
P2nucs = []
for  n, eos in enumerate (eos_list):
    # This is too silly
    P2nucs.append(np.log10(lalsim.SimNeutronStarEOSEnergyDensityOfPressure(2*nuc_rho_si, eos.eos)/10./ccgs**2)) # In cgs
    if max_masses[n] > 1.405:
        R1p4s.append(lalsim.SimNeutronStarRadius(M1p4s, eos.family)/1000.)
    else:
         print("failed to produce a reasonable max mass")
#Loves = np.array([lalsim.SimNeutronStarLoveNumberK2(M1p4s, eos.family) for eos in eos_list])

# Mass plot
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


# Pressure(2 rho_nuc)
plt.figure()
plt.hist(P2nucs, bins=30)
plt.xlabel(r"$\log_{10}(P(2\rho_{nuc})/c^2$)")
plt.ylabel("frequency")
plt.savefig("better_prior_P2nuc.png")
plt.show()

