#!/usr/bin/env python
'''
Generates a figure ('shock_profiles.png') that profiles various key variables. The temperature,
pressure, and radiation energy density should be proportional to exp(-sqrt(tau)) where tau is
the optical depth measured from the shock front. This is true at least in the case of the
subcritical shock. Unsure about supercritical shocks.

This script can be run without any command-line arguments, but requires that yt version 3.0+
be installed.
'''
import numpy as np
import matplotlib.pyplot as plt
from yt import *  # Requires yt version 3.0+
import sys
import glob

# Font size
from matplotlib import rcParams as rcp
rcp['font.size'] = 12
rcp['xtick.labelsize'] = 12
rcp['ytick.labelsize'] = 12

pltfiles = np.sort(glob.glob('subcritical_shock_hdf5_plt_cnt_*'))

fig = plt.figure()
ax1 = plt.subplot(221)
colors = palette = [\
               '#EAD846',\
               '#6F0989',\
               '#D9712A',\
               '#97C9E4',\
               '#B82035',\
               '#C3C385',\
               '#62AC49',\
               '#CE81AD',\
               '#476CB3',\
               ]

# Take the last plotfile and profile it
pf = load(pltfiles[-1])
domain = pf.domain_width[0].value

A = np.array([0.0,0,0])*domain
B = np.array([1.0,0,0])*domain

ray = pf.ray(A,B)
v   = pf.parameters['sim_velx']
t   = pf.current_time.value

x    = ray['x'].value
dx   = ray['dx'].value
xx   = x - v*t
tele = ray['tele'].value               # Electron temperature
tion = ray['tion'].value               # Ion temperature. Should have T_ele = T_ion
trad = ray['trad'].value               # Radiation temperature
temp = ray['temp'].value               # Value of 'temp' variable. Depends on runtime parameter settings
erad = ray['erad'].value               # Specific radiation energy density
try:
    kapp = ray['kapp'].value           # Opacity in cm^2/g
except:
    kapp = pf.parameters['op_absorbconst']
dens = ray['dens'].value               # Mass density
pres = ray['pres'].value               # Hydrostatic pressure
tau  = np.cumsum(kapp*dx)              # Optical depth measured from x = 0 in lab frame

shock_ind   = tele.argmax()            # Array index of shock position
xshock      = x[shock_ind]             # Shock coordinates
shock_speed = x[shock_ind]/t           # Speed of shock in lab frame
u           = xx[shock_ind]/t          # Speed of incoming material / gas velocity in shock frame
taus        = tau - tau[shock_ind]     # Optical depth as measured from the shock
taus[taus < 0.0] = 0.0
print 'Shock is traveling at {0} km/s to the right'.format(shock_speed/1.0e5)

plt.semilogy(x,tele,color=colors[5])
plt.semilogy(x,tion,color=colors[8])
plt.semilogy(x,trad,color=colors[6])
fit = 0.8*tele[shock_ind+10]*np.exp(-np.sqrt(3)*taus)
plt.semilogy(x[taus>0],fit[taus>0],'--',color=colors[4])

plt.xlabel(r'Distance [cm]')
plt.ylabel(r'Temperature [K]')

ax1.set_xlim([0, 0.5*domain])
ax1.set_ylim([8,1000])
plt.legend(['ele','ion','rad',r'e$^{-\sqrt{3}\tau}$'])

ax2 = plt.subplot(222)
fit = 0.8*pres[shock_ind+10]*np.exp(-np.sqrt(3)*taus)
plt.semilogy(x,pres)
plt.semilogy(x[taus>0],fit[taus>0],'--')

plt.xlabel(r'Distance [cm]')
plt.ylabel(r'Pressure [dyn/cm$^2$]')

ax2.set_xlim([0, 0.5*domain])
ax2.set_ylim([1,1000])

ax3 = plt.subplot(223)
plt.semilogy(x,dens)
fit = dens[shock_ind]*np.exp(-np.sqrt(3)*taus)
plt.semilogy(x[taus>0],fit[taus>0],'--')

plt.xlabel(r'Distance [cm]')
plt.ylabel(r'Density [g/cm$^3$]')

ax3.set_xlim([0, 0.5*domain])
ax3.set_ylim([7e-10,1e-8])

ax4 = plt.subplot(224)
plt.semilogy(x,erad)
fit = 10.0*erad[shock_ind]*np.exp(-np.sqrt(3)*taus)
plt.semilogy(x[taus>0],fit[taus>0],'--')

plt.xlabel(r'Distance [cm]')
plt.ylabel(r'Radiation Density [erg/g]')

ax4.set_xlim([0, 0.5*domain])


# Adjust the labels
axs = [ax1, ax2, ax3, ax4]
for ax in axs:
    ax.tick_params(axis='x', pad = 10)
    ax.tick_params(axis='y', pad = 10)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10

plt.tight_layout()
plt.savefig('shock_profiles.png')
