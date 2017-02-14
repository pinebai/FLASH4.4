#if 0
This file contains the indices of the various variables used in 
the EOS routines.  The first seven indices (EOS_VARS) are basic properties
whereas the last eleven indices (EOS_VARS+1 to EOS_NUM) represent derived 
or optionally calculated quantities.

The derived variables are used in the mask required by Eos.F90, and
values for them are returned by Eos in the array eosData when the
corresponding mask element is true.

NOTE1 -- originally both EOS_ENER and EOS_EINT were included here.  
EOS_ENER was the variable (mostly) used internally, and EOS_EINT was basically 
  ignored.  This was very confusing, as the EOS routines work with internal
  entergy, which corresponds to EINT_VAR in the unk array.  Hence, on
7/20/2006 LBR removed all references to EOS_ENER by changing them to EOS_EINT.

NOTE2 -- the equivalent Flash2 variable is given in [brackets] after the 
  description below

quantity EOS_VARS 9
quantity EOS_DERIVS 13
quantity EOS_NUM 22

NOTE3a -- kinetic energy EKIN is not actually used by any Eos implementation,
but a place for it reserved among the "basic" variables so that an Eos-wrapping
function can conveniently store the kinetic energy derived from grid
datastructure variables across Eos() calls.

basic EOS_PRES  Pressure [p]
basic EOS_DENS  Density  [rho]
basic EOS_TEMP  Temperature [temp]
basic EOS_GAMC  Adiabatic index, the Gamma1 of Chandrasekhar [gamc]
basic EOS_EINT  Total internal energy (equivalent to EINT_VAR) [ei]
basic EOS_ABAR  Average mass of the nuclei
basic EOS_ZBAR  Average proton number
basic EOS_ENTR  Entropy 
basic EOS_EKIN  used to hold kinetic energy for internal Eos use
derived EOS_DPT Derivative of pressure wrt temperature
derived EOS_DPD Derivative of pressure wrt density
derived EOS_DET Derivative of internal energy wrt temperature
derived EOS_DEA Derivative of internal energy wrt atomic mass
derived EOS_DEZ Derivative of internal energy wrt atomic number/charge
derived EOS_DED Derivative of internal energy wrt density
derived EOS_DST Derivative of entropy wrt temperature  
derived EOS_DSD Derivative of entropy wrt density      
derived EOS_CV  Specific heat at constant volume
derived EOS_CP  Specific heat at constant pressure
derived EOS_PEL Electron pressure
derived EOS_NE  Electron number density
derived EOS_ETA Electron degeneracy parameter (chemical potential / k_b*T)

#endif

#define EOS_BEGIN 1
#define EOS_VARS 12
#define EOS_DERIVS 23
#define EOS_NUM 35
#define EOS_PRES 1
#define EOS_DENS 2
#define EOS_EINT 3
#define EOS_TEMP 4
#define EOS_GAMC 5
#define EOS_ABAR 6
#define EOS_ZBAR 7
#define EOS_ENTR 8
#define EOS_EKIN 9

#define EOS_TEMPION 10
#define EOS_TEMPELE 11
#define EOS_TEMPRAD 12

#define EOS_DPT 13
#define EOS_DPD 14
#define EOS_DET 15
#define EOS_DED 16
#define EOS_DEA 17
#define EOS_DEZ 18
#define EOS_DST 19
#define EOS_DSD 20
#define EOS_CV  21
#define EOS_CP  22
#define EOS_PEL 23
#define EOS_NE  24
#define EOS_ETA 25

#define EOS_EINTION 26
#define EOS_EINTELE 27
#define EOS_EINTRAD 28
#define EOS_PRESION 29
#define EOS_PRESELE 30
#define EOS_PRESRAD 31

#define EOS_ENTRELE 32
#define EOS_ENTRRAD 33

#define EOS_CVION   34
#define EOS_CVELE   35

#define N_EOS_TEMP 3

#if 0
 This section defines the constants to identify different equations
 of states that are included in the distribution.
 The "USER" types can be used for user-defined specialized implementations.
#endif

#define EOS_GAM 321
#define EOS_HLM 452
#define EOS_MGAM 645
#define EOS_MTMP 325
#define EOS_TAB 624
#define EOS_MTYPE 625

#define EOS_USERTYPE1 1001
#define EOS_USERTYPE2 1002


#define EOS_TABULAR_Z 1
#define EOS_TABULAR_E 2 
#define EOS_TABULAR_C 3
#define EOS_TABULAR_P 4
#define EOS_TABULAR_S 5
#define EOS_APPROX_KIN 40


#define EOS_LOGLEVEL_WARN_ANY      100
#define EOS_LOGLEVEL_WARN_DATA     300
#define EOS_LOGLEVEL_WARN_ALLPROCS 430
#define EOS_LOGLEVEL_INFO_DATA     450
#define EOS_LOGLEVEL_INFO_ALLPROCS 900
#define EOS_LOGLEVEL_INFO_ALL     1000
