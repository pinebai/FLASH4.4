#!/usr/bin/env python

from math import *
from consts import *

N    = 10000
zmin = 0.0*cmpc
zmax = 1000.*cmpc
dz   = (zmax - zmin) / N

T = 1e4
mu = 0.609
gamma = 1.0001
cs2 = gamma*kB*T/mu/mH
rho0 = 1*mH

def sech(x):
    return 2.*exp(x)/(1.+exp(2.*x))

def cosh(x):
    return (1.+exp(2.*x))/(2.*exp(x))

print N

z = zmin
while z < zmax:
    rho = rho0*sech(sqrt(2*pi*G*rho0/cs2)*z)**2
    p = rho*cs2/gamma
    v = 0.0
    x = 1.0
    Phi = 2*cs2*log(cosh(sqrt(2*pi*G*rho0/cs2)*z))
    print '% 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e ' % (z, rho, p, v, x, Phi)

    z += dz


