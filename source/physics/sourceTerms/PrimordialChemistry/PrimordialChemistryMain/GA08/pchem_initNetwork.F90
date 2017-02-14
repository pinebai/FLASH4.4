!!****ih* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_initNetwork
!!
!! NAME
!!
!! pchem_initNetwork
!!
!!
!! SYNOSIS
!!
!! pchem_initNetwork()
!!
!!
!! DESCRIPTION
!!
!! Initializes variables and integers used in this network.
!!
!! Called from PrimordialChemistry_init 
!!
!! this is going to include all the reactions from S.C.O Glover and T. Abel 
!! well, most of them
!!
!! ions: H+,H,H-,D+,D,D-,He++,He+,HE,H2+,H2,HD+,HD,e-
!!
!!***

subroutine pchem_initNetwork

    use PrimordialChemistry_data
    use pchem_data 

#include "Flash.h"
  implicit none

  integer :: i

  !..Zero all the isotope pointers
  do i=i,nisotp
     isotp(i) = 0
  enddo

  !..zero the steps taken
  xoktot = 0.0e0
  xbadtot = 0.0e0
  xkburn = 0.0e0

  !..set the id numbers of the ions used
  iH    = 1
  iHP   = 2
  iHM   = 3
  iD    = 4
  iDP   = 5
  iDM   = 6
  iHE   = 7
  iHEPP = 8
  iHEP  = 9
  iH2P  = 10
  iH2   = 11
  iHDP  = 12
  iHD   = 13
  iD2   = 14
  iD2P  = 15
  iELEC = 16

  
  ionam(iHP)   = 'HP'
  ionam(iH)    = 'H'
  ionam(iHM)   = 'HM'
  ionam(iDP)   = 'DP'
  ionam(iD)    = 'D'
  ionam(iDM)   = 'DM'
  ionam(iHEPP) = 'HEPP'
  ionam(iHEP)  = 'HEP'
  ionam(iHE)   = 'HE'
  ionam(iH2P)  = 'H2P'
  ionam(iH2)   = 'H2'
  ionam(iHDP)  = 'HDP'
  ionam(iHD)   = 'HD'
  ionam(iD2)   = 'D2'
  ionam(iD2P)  = 'D2P'
  ionam(iELEC) = 'ELEC'

  !..set the number of nucleons in each element
  aion(iH)     = 1.0e0
  aion(iHP)      = 0.999451
  aion(iHM)     = 1.000549
  aion(iD)     = 2.0e0
  aion(iDP)      = 1.999451
  aion(iDM)     = 2.000549
  aion(iHE)   = 4.0e0
  aion(iHEPP)    = 3.998902
  aion(iHEP)     = 3.999451
  aion(iH2P)    = 1.999451
  aion(iH2)     = 2.0e0
  aion(iHDP)    = 2.999451
  aion(iHD)     = 3.0e0
  aion(iELEC)   = 0.000549
  aion(iD2)     = 4.000
  aion(iD2P)    = 3.999451

  do i=1,NSPECIES
     aioninv(i) = 1.0e0/aion(i)
  enddo

  !..set the number of protons in each element
  zion(iHP)     = 1.0e0
  zion(iH)      = 1.0e0
  zion(iHM)     = 1.0e0
  zion(iDP)     = 1.0e0
  zion(iD)      = 1.0e0
  zion(iDM)     = 1.0e0
  zion(iHEPP)   = 2.0e0
  zion(iHEP)    = 2.0e0
  zion(iHE)     = 2.0e0
  zion(iH2P)    = 2.0e0
  zion(iH2)     = 2.0e0
  zion(iHDP)    = 2.0e0
  zion(iHD)     = 2.0e0
  zion(iELEC)   = 0.0e0
  zion(iD2)     = 2.0e0
  zion(iD2P)    = 2.0e0
  do i=1,NSPECIES
     zionsq(i) = zion(i) * zion(i)
  enddo

  !..set the binding energy of each element
  !..I will set this as zero for now
  !..Hmm, Binding energies seem to be in MeV
  bion(iH)    = 00.000000e-6
  bion(iHP)   = -13.598440e-6
  bion(iHM)   = 0.750000e-6
  bion(iD)    = 00.000000e-6
  bion(iDP)   = -13.598440e-6
  bion(iDM)   = 0.750000e-6
  bion(iHEPP) = -79.005100e-6
  bion(iHEP)  = -24.587410e-6
  bion(iHE)   = 00.000000e-6
  bion(iH2P)  = -10.947100e-6
  bion(iH2)   = 4.477880e-6
  bion(iHDP)  = -10.947100e-6
  bion(iHD)   = 4.477880e-6
  bion(iELEC) = 00.000000e-6
  bion(iD2)   = 4.477880e-6
  bion(iD2P)  = -10.947100e-6

  !..set the id numbers of the reaction rates
  !..I am going to create a shorthand here and explain
  !..each one. Basically the idea will be each one will start
  !.. with 'ir' first, followed by the reactants and then the 
  !..products, but I will not put down electrons or gammas

  iR001 = 1
iR002 = 2
iR003 = 3
iR004 = 4
iR005 = 5
iR006 = 6
iR007 = 7
iR008 = 8
iR009 = 9
iR010 = 10
iR011 = 11
iR012 = 12
iR013 = 13
iR014 = 14
iR015 = 15
iR016 = 16
iR017 = 17
iR018 = 18
iR019 = 19
iR020 = 20
iR021 = 21
iR022 = 22
iR023 = 23
iR024 = 24
iR025 = 25
iR026 = 26
iR027 = 27
iR028 = 28
iR029 = 29
iR030 = 30
iR031 = 31
iR032 = 32
iR033 = 33
iR034 = 34
iR035 = 35
iR036 = 36
iR037 = 37
iR038 = 38
iR039 = 39
iR040 = 40
iR041 = 41
iR042 = 42
iR043 = 43
iR044 = 44
iR045 = 45
iR046 = 46
iR047 = 47
iR048 = 48
iR049 = 49
iR050 = 50
iR051 = 51
iR052 = 52
iR053 = 53
iR054 = 54
iR055 = 55
iR056 = 56
iR057 = 57
iR058 = 58
iR059 = 59
iR060 = 60
iR061 = 61
iR062 = 62
iR063 = 63
iR064 = 64
iR065 = 65
iR066 = 66
iR067 = 67
iR068 = 68
iR069 = 69
iR070 = 70
iR071 = 71
iR072 = 72
iR073 = 73
iR074 = 74
iR075 = 75
iR076 = 76
iR077 = 77
iR078 = 78
iR079 = 79
iR080 = 80
iR081 = 81
iR082 = 82
iR083 = 83
iR084 = 84
iR085 = 85
iR086 = 86
iR087 = 87
iR088 = 88
iR089 = 89
iR090 = 90
iR091 = 91
iR092 = 92
iR093 = 93
iR094 = 94
iR095 = 95
iR096 = 96
iR097 = 97
iR098 = 98
iR099 = 99
iR100 = 100
iR101 = 101
iR102 = 102
iR103 = 103
iR104 = 104
iR105 = 105
iR106 = 106
iR107 = 107
iR108 = 108
iR109 = 109
iR110 = 110
iR111 = 111
iR112 = 112
iR113 = 113
iR114 = 114
iR115 = 115
iR116 = 116
iR117 = 117
iR118 = 118
iR119 = 119
iR120 = 120
iR121 = 121
iR122 = 122
iR123 = 123
iR124 = 124


 return
end subroutine pchem_initNetwork
