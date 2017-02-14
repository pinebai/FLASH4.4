#!/usr/bin/env python

from math import *

# based on Bonnor, W. B., 1956, MNRAS, 116, 351

# physical constants and unit coversions
mH     = 1.6605387e-24
kB     = 1.3806503e-16
G      = 6.673e-08
MSun   = 1.9891e+33
cmpc   = 3.0856776e+18

# BES parameters
mfac   = 1.5    # mass factor - to scale BES to a higher mass to get the collapse
M      = MSun   # mass of the BES
T      = 10.0   # temperature in the BES
xi0    = 6.0    # dimensionless radius at the BES edge
T_amb  = 1.0e4  # temperature in the ambient medium
m      = 2*mH   # mass of a particle
N      = M/m    # number of particles
beta   = kB*T / (4*pi*G*m) # see B56

Nstep  = 1e4    # number of integration steps
Nprint = 1e3    # number of printed lines
rend   = 0.2*cmpc  # maximum radius


def p(T, m, N, xi, psi, dpsidxi):
    return (kB*T)**4/(4*pi*G**3*m**6*N**2) * xi**4 * dpsidxi**2 * exp(-psi)

def V(T, m, N, xi, psi, dpsidxi):
    return 4*pi/3 * (G*m*m*N/(kB*T*xi*dpsidxi+1e-99))**3

def equ31(xi, psi, dpsidxi):
    return 1.0 - 0.5*exp(psi) * dpsidxi**2

def equ32(xi, psi, dpsidxi):
    return 1.0 - exp(psi)/(xi+1e-99) * dpsidxi


def rk4(xi, psi, eta, dxi):
    k1_psi = dxi * eta
    k1_eta = dxi * (exp(-psi) - 2*eta/(xi+1e-99))

    k2_psi = dxi * (eta + 0.5*k1_eta)
    k2_eta = dxi * (exp(-(psi + 0.5*k1_psi)) - 2*(eta + 0.5*k1_eta)/(xi + 0.5*dxi))

    k3_psi = dxi * (eta + 0.5*k2_eta)
    k3_eta = dxi * (exp(-(psi + 0.5*k2_psi)) - 2*(eta + 0.5*k2_eta)/(xi + 0.5*dxi))

    k4_psi = dxi * (eta + k3_eta)
    k4_eta = dxi * (exp(-(psi + k3_psi)) - 2*(eta + k3_eta)/(xi + dxi))

    psi_new = psi + (k1_psi/6.0 + k2_psi/3.0 + k3_psi/3.0 + k4_psi/6.0)
    eta_new = eta + (k1_eta/6.0 + k2_eta/3.0 + k3_eta/3.0 + k4_eta/6.0)
    return (psi_new, eta_new)



# find central density
xistart = 0.0
xiend   = xi0
dxi     = (xiend-xistart)/Nstep
psi     = 0.0
dpsidxi = 0.0
xi      = xistart
while xi < xi0:
    (psi, dpsidxi) = rk4(xi, psi, dpsidxi, dxi)
    xi += dxi
rho0     = (4*pi)**2/(m*N)**2 * beta**3 * xi**4 * dpsidxi**2

# conditions at the BES edge
r_edge   = sqrt(beta/rho0)*xi
rho_edge = rho0*exp(-psi)
p_edge   = rho_edge*kB*T/m
Phi_edge = -kB*T/m * log(rho_edge)
Phi0     = kB*T/m * log(rho_edge) - G*M/r_edge


# initial conditions for the integration
xistart = 0.0
xiend   = xi0
dxi     = (xiend-xistart)/Nstep

psi     = 0.0
dpsidxi = 0.0
xi      = xistart

dprint  = rend/Nprint
rprint  = 0.0

r   = 0.0
rho = rho0
p   = rho*kB*T/m
v   = 0.0
x   = 1.0
Phi = -kB*T/m * log(rho) + Phi0

# solution arrays
Nr = 0
(rarr, rhoarr, parr, varr, xarr, Phiarr) = ([], [], [], [], [], [])

# add the first point
rarr.append(r)
rhoarr.append(rho)
parr.append(p)
varr.append(v)
xarr.append(x)
Phiarr.append(Phi)
Nr += 1

while r < rend:
    (psi, dpsidxi) = rk4(xi, psi, dpsidxi, dxi)
    xi += dxi
    r   = sqrt(beta/rho0)*xi
    if r < r_edge: 
        rho = rho0*exp(-psi)
        p   = rho*kB*T/m
        x = 1.0
        Phi = -kB*T/m * log(rho) + Phi0

    else:
        p = p_edge
        rho = p_edge*m/(kB*T_amb)
        x = 0.0
        Phi = -G*M/r


    if (r - rprint) >= dprint:
        rarr.append(r)
        rhoarr.append(rho)
        parr.append(p)
        varr.append(v)
        xarr.append(x)
        Phiarr.append(Phi)
        Nr += 1
        rprint = r

# print solution
print Nr
for i in range(Nr):
    print "% 15.8e  % 15.8e  % 15.8e  % 15.8e  % 15.8e  % 15.8e" \
    % (rarr[i], rhoarr[i]*mfac, parr[i]*mfac, varr[i], xarr[i], Phiarr[i]*mfac)

print
print "# r               rho              p                v                x                Phi"
print "# M        = %12.5e MSun = %12.5e g" % (M/MSun*mfac, M*mfac)
print "# rho0     = %12.5e g/cm^3 = %12.5e cm^-3" % (rho0*mfac, rho0/m*mfac)
print "# r_edge   = ", r_edge, "cm = ", r_edge/cmpc, " pc"
print "# rho_edge = %12.5e g/cm^3 = %12.5e cm^-3" % (rho_edge*mfac, rho_edge/m*mfac)
print "# p_edge   = ", p_edge*mfac, " dyne/cm^2"
print "# Phi_edge = ", Phi_edge*mfac, " cm^2/s^2"


