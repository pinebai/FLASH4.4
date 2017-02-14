#if 0
This file contains symbols for integer indices that represent the various ROLES
in which variables in Grid data structures, including UNK, GRIDVAR, and face variables,
can be used.  These symbols therefore provide a generalization of the set of UNK
variables that have traditionally been used as input to, and output from, Eos
calls via Eos_wrapped.

Some of the corresponding variables in Grid data structures get directly
copied to corresponding basic Eos data structure components passed in calls to Eos(),
and/or copied back from corresponding basic Eos data structure components after
Eos calls; see Eos.h for the list of defined Eos data structure components.
For other variables in Grid data structures that correspond to Eos_map roles,
the relationship to their related Eos data structure components is less direct.
For example, the Grid data structures that correspond to Eos_map roles EOSMAP_VEL[XYZ]
are used as inputs to compute kinetic energy for the purpose of determining the value
of the EOS_EINT basic Eos data structure component when calling Eos().

See Eos_getData.F90 and Eos_putData.F90 for how these Eos_map roles are actually used.

EOSMAP_PRES  Pressure    (normally PRES_VAR)
EOSMAP_DENS  Density     (normally DENS_VAR)
EOSMAP_TEMP  Temperature (normally TEMP_VAR)
EOSMAP_GAMC  Adiabatic index, Chandrasekhars Gamma1 (normally GAMC_VAR)
EOSMAP_GAME  Adiabatic index, another of Chandrasekhars Gammas (normally GAME_VAR)
EOSMAP_EINT  Total internal energy (normally EINT_VAR)
EOSMAP_ENER  Total energy, internal plus kinetic (normally ENER_VAR)
EOSMAP_SUME  Used in Ye-based Eos to derive average mass of the nuclei
EOSMAP_YE    Used in Ye-based Eos to derive Average proton number
EOSMAP_ENTR  Entropy, normally output only
EOSMAP_VELX  Velocity component, used to compute kinetic energy, input only
EOSMAP_VELY  Velocity component, used to compute kinetic energy, input only
EOSMAP_VELZ  Velocity component, used to compute kinetic energy, input only

This can easily be added if the need arises:
EOSMAP_EKIN  used to hold kinetic energy for internal Eos use
#endif


#define EOSMAP_BEGIN 1
#define EOSMAP_NUM_ROLES 27

#define EOSMAP_PRES 1
#define EOSMAP_DENS 2
#define EOSMAP_EINT 3
#define EOSMAP_TEMP 4
#define EOSMAP_GAMC 5
#define EOSMAP_GAME 6
#define EOSMAP_ENER 7
#define EOSMAP_VELX 8
#define EOSMAP_VELY 9


#define EOSMAP_VELZ 10
#define EOSMAP_ENTR 11
#define EOSMAP_SUMY 12
#define EOSMAP_YE   13

#define EOSMAP_PRES1 14
#define EOSMAP_PRES2 15
#define EOSMAP_PRES3 16
#define EOSMAP_EINT1 17
#define EOSMAP_EINT2 18
#define EOSMAP_EINT3 19
#define EOSMAP_TEMP1 20
#define EOSMAP_TEMP2 21
#define EOSMAP_TEMP3 22
#define EOSMAP_E1 23
#define EOSMAP_E2 24
#define EOSMAP_E3 25

#define EOSMAP_SELE 26
#define EOSMAP_SRAD 27

#define EOS_IN 1
#define EOS_OUT 2
