import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['figure.figsize'] = (9.7082039325, 6.0)
matplotlib.rcParams['xtick.labelsize'] = 20.0
matplotlib.rcParams['ytick.labelsize'] = 20.0
matplotlib.rcParams['axes.labelsize'] = 25.0
matplotlib.rcParams['legend.fontsize'] = 16.0
matplotlib.rcParams['font.family']= 'Times New Roman'
matplotlib.rcParams['font.sans-serif']= ['Bitstream Vera Sans']
matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['mathtext.fontset']= 'stixsans'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True

import lal
import scipy
import random
from random import randint
import csv
from scipy.interpolate import interp1d

c_cgs=lal.C_SI*100
rhonuc=2.8e14


fig = plt.figure(figsize=(9.7082039325, 6))

path="../../prior_quantiles_cs2.csv"
columns = open(path, 'r').readline().strip().split(',')
data = np.loadtxt(path, delimiter=',', skiprows=1)
baryon_density = [float(_.split('=')[1][:-1]) for _ in columns[1:]]
cs2s50 = data[50,:]
cs2s5 = data[5,:]
cs2s95 = data[95,:]

plt.fill_between(baryon_density,cs2s5[1:],cs2s95[1:],color='c',alpha=0.05)
plt.plot(baryon_density,cs2s5[1:],c='k',lw=2.5,label='Prior')
plt.plot(baryon_density,cs2s95[1:],c='k',lw=2.5)

path="../../no_j0740_quantiles_cs2.csv"
columns = open(path, 'r').readline().strip().split(',')
data = np.loadtxt(path, delimiter=',', skiprows=1)
baryon_density = [float(_.split('=')[1][:-1]) for _ in columns[1:]]
cs2s50 = data[50,:]
cs2s5 = data[5,:]
cs2s95 = data[95,:]

plt.fill_between(baryon_density,cs2s5[1:],cs2s95[1:],color='c',alpha=0.2)
plt.plot(baryon_density,cs2s5[1:],c='c',lw=2.5,label='PSR + GW + J0030')
plt.plot(baryon_density,cs2s95[1:],c='c',lw=2.5)

path="../../all_miller_quantiles_cs2.csv"
columns = open(path, 'r').readline().strip().split(',')
data = np.loadtxt(path, delimiter=',', skiprows=1)
baryon_density = [float(_.split('=')[1][:-1]) for _ in columns[1:]]
cs2s50 = data[50,:]
cs2s5 = data[5,:]
cs2s95 = data[95,:]

plt.fill_between(baryon_density,cs2s5[1:],cs2s95[1:],alpha=0.5)
#plt.plot(baryon_density,cs2s50[1:],c='b')
plt.plot(baryon_density,cs2s5[1:],c='b',lw=2.5,label='PSR + GW + J0030 + J0740')
plt.plot(baryon_density,cs2s95[1:],c='b',lw=2.5)

plt.axvline(rhonuc,c='k')
plt.axvline(2*rhonuc,c='k')
plt.axvline(6*rhonuc,c='k')
plt.axhline(1/3, c='k')

plt.tick_params(direction='in')

plt.xscale('log')

plt.tight_layout()
plt.grid(alpha=0.5)
plt.ylabel('$ c_s^2/c^2$')
plt.xlabel('$\\rho \,(\mathrm{g/cm}^3)$')
plt.xlim(4e13,3e15)
plt.ylim(0.0,1.0)
plt.legend(frameon=True,fancybox=True,framealpha=1,fontsize=18)
plt.text(3e14,0.1,'$\\rho_{\mathrm{nuc}}$',fontsize=28,rotation=90)
plt.text(6e14,0.1,'$2\\rho_{\mathrm{nuc}}$',fontsize=28,rotation=90)
plt.text(18e14,0.1,'$6\\rho_{\mathrm{nuc}}$',fontsize=28,rotation=90)
plt.text(5e13,0.37,'$c_s^2 = 1/3$',fontsize=22)

plt.show()
fig.savefig("../../plot_cs2-rho_quantiles_comparison.pdf",bbox_inches='tight')
