#!/usr/bin/env python
'''
Draws the profile of the matter and radiation temperature through the shock and compares
the values to some semi-analytic estimates from Mihalas & Mihalas (1984). See also Klassen et
al. (2014).

The agreement between measured values and the semi-analytic estimates is not very good.
The post-shock gas temperature is too cold. The temperature in the spike is too high. The
temperature in the radiative precursor, just before the spike, is also higher than expected.

We believe these inconsistences are the result of coupling between the matter species, and
between the radiation and matter.

This script should run without needing any command-line arguments. The main requirement is
that yt version 3.0+ be installed on the system.

There is one boolean called "subcritical" that must be set to True or False, depending on
whether the simulation was for the subcritical (True) or supercritical (False) shock.
'''
import numpy as np
import matplotlib.pyplot as plt
from yt import *   # requires yt version 3.0+
import sys
import glob

# Font size
from matplotlib import rcParams as rcp
rcp['font.size'] = 16
rcp['xtick.labelsize'] = 20
rcp['ytick.labelsize'] = 20

subcritical = True

if subcritical:
    pltfiles = np.sort(glob.glob('subcritical_shock_hdf5_plt_cnt_*'))
else:
    pltfiles = np.sort(glob.glob('supercritical_shock_hdf5_plt_cnt_*'))

fig = plt.figure()
ax = plt.subplot(111)
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

for i in range(len(pltfiles)):
    pf = load(pltfiles[i])
    domain = pf.domain_width[0].value

    A = np.array([0.01,0,0])*domain
    B = np.array([0.5,0,0])*domain

    ray = pf.ray(A,B)
    v   = pf.parameters['sim_velx']
    t   = pf.current_time.value

    x    = ray['x'].value
    xx   = x - v*t
    tele = ray['tele'].value
    tion = ray['tion'].value
    trad = ray['trad'].value
    temp = ray['temp'].value

    plt.plot(x,tele,color=colors[4])
    shock_ind = tele.argmax()
    shock_speed = x[shock_ind]/t
    u = xx[shock_ind]/t
    print 'Shock is traveling at {0} km/s to the right'.format(shock_speed/1.0e5)
    plt.plot(x,tion,color=colors[8])
    plt.plot(x,trad,color=colors[6])

    plt.xlabel('Distance [cm]')
    plt.ylabel('Temperature [K]')

ax.set_xlim([0, 0.5*domain])
plt.legend(['ele','ion','rad'])

# Adjust the labels
ax.tick_params(axis='x', pad = 10)
ax.tick_params(axis='y', pad = 10)
ax.xaxis.labelpad = 10
ax.yaxis.labelpad = 10

plt.legend(['Electrons','Ions','Radiation'],loc=1)
plt.savefig('shock_temperatures.png')

dx = pf.index.get_smallest_dx().in_units('cm')

# Make a figure of just the last plotfile
plt.clf()

fig = plt.figure()
ax = plt.subplot(111)
pf = load(pltfiles[i])

# Extract the useful parameters and constants
f = open('shock_simulation.info','w')
rho = 7.78e-10
kb = 1.38065E-16
mH = 1.67262E-24
sigma = 5.67040E-05
gamma = pf.h.parameters['gamma']
try:
    mu = pf.parameters['eos_mu_mol']
except:
    print 'eos_mu_mol not set. Using eos_singleSpeciesA'
    mu = pf.parameters['eos_singlespeciesa']
# Gas constant
R = 8.31447E+07
mu = 0.5
R = kb/mu/mH 
print 'mu_mol = ', mu
print 'R = ',R
vars = {'rho':rho, 'R':R, 'gamma':gamma,'mu':mu}
for var in vars:
    f.write('\n '+var+' = '+str(vars[var]))

ray = pf.h.ray(A,B)
v   = pf.parameters['sim_velx']
t   = pf.current_time.value

x    = ray['x'].value
xx   = x - v*t
tele = ray['tele'].value
tion = ray['tion'].value
trad = ray['trad'].value
temp = ray['temp'].value
tele = ray['tele'].value
shock_ind = tele.argmax()
u = xx[shock_ind]/t

plt.plot(x,temp,'o',color=[0.4,0.4,0.4])
plt.plot(x,tele)
plt.plot(x,tion)
p1, = plt.plot(x,temp,color=colors[4])
p2, = plt.plot(x,trad,'--',color=colors[6])

plt.xlabel('Distance [cm]')
plt.ylabel('Temperature [K]')

ax.set_xlim([0, B[0]])

# Adjust the labels
ax.tick_params(axis='x', pad = 10)
ax.tick_params(axis='y', pad = 10)
ax.xaxis.labelpad = 10
ax.yaxis.labelpad = 10

plt.legend([p1, p2], ['Matter','Radiation'],loc=1)
if subcritical:
    plt.savefig('subcritical_shock.png')
else:
    plt.savefig('supercritical_shock.png')

# Compare to analytic values
T2 = 2.0*(gamma - 1.0)*u**2 / (R*(gamma+1.0)**2)
T2sim = temp[0]
Tm = (gamma - 1.0)/(rho*abs(u)*R) * (2.0*sigma*T2sim**4)/(np.sqrt(3))
Tp = T2sim + (3.0 - gamma)/(gamma + 1.0) * Tm

if subcritical:
    print 'The post-shock gas temperature is supposed to be:'
    print T2
    print 'Whereas we measure:'
    print temp[0]
    f.write('\n\nThe post-shock gas temperature is supposed to be: {0}'.format(T2))
    f.write('\nWhereas we measure: {0}'.format(temp[0]))
    print 'The preshock gas temperature in the radiative precursor should be:'
    print Tm
#    print 'Whereas we measure:'
#    print temp[np.argmax(temp):np.argmax(temp)+10]
#
#    Tmsim = input('Please enter the preshock temp read from the array above:  ')
    Tp = T2sim + (3.0 - gamma)/(gamma + 1.0) * Tm
    #Tp = T2sim + (3.0 - gamma)/(gamma + 1.0) * Tmsim
    print 'The temperature spike should be:'
    print Tp
    print 'Whereas we measure:'
    print max(temp)
else:
    # Supercritical case
    Tp = (3.0 - gamma)*T2sim
    print 'For the supercritical case:'
    print 'The temperature spike should be:'
    print Tp
    print 'Whereas we measure:'
    print max(tion)

print '\nOther information:'

print 'Smallest dx = {0}'.format(dx)
print 'Time = {0}'.format(t)
f.write('\n dx   = {0} '.format(dx))
f.write('\n time = {0} '.format(t))
f.write('')
f.close()
