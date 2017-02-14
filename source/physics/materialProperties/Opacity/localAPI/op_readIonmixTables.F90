!!****if* source/physics/materialProperties/Opacity/localAPI/op_readIonmixTables
!!
!! NAME
!!
!!  op_readIonmixTable
!!
!! SYNOPSIS
!!
!!  call op_readIonmixTables (character (in) :: tableName (len=80),
!!                            logical   (in) :: needPATable,
!!                            logical   (in) :: needPETable,
!!                            logical   (in) :: needROTable,
!!                            integer   (in) :: indexPA,
!!                            integer   (in) :: indexPE,
!!                            integer   (in) :: indexRO)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an IONMIX datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays:
!!
!!             op_PlanckAbsorptionTables (t,d,g,indexPA)  (in cm^2/g)
!!             op_PlanckEmissionTables   (t,d,g,indexPE)  (in cm^2/g)
!!             op_RosselandTables        (t,d,g,indexRO)  (in cm^2/g)
!!
!!    where:   t = temperature (eV) index
!!             d = ion number density (# ions/cm^3) index
!!             g = energy group (eV) index
!!       indexXX = table counting index for Opacity kind XX (XX = PA,PE,RO)
!!
!!  The number of temperature and density indices will be stored in the 1-dimensional arrays: 
!!
!!             op_nstepsDensityPA  (indexPA)
!!             op_nstepsDensityPE  (indexPE)
!!             op_nstepsDensityRO  (indexRO)
!!
!!             op_nstepsTemperaturePA  (indexPA)
!!             op_nstepsTemperaturePE  (indexPE)
!!             op_nstepsTemperatureRO  (indexRO)
!!
!!  The actual temperatures and densities will be stored in the 2-dimensional arrays:
!!
!!           op_tableDensityPA     (d,indexPA) ; d = 1,op_nstepsDensityPA (indexPA)       (in # ions/cm^3)
!!           op_tableDensityPE     (d,indexPE) ; d = 1,op_nstepsDensityPE (indexPE)       (in # ions/cm^3)
!!           op_tableDensityRO     (d,indexRO) ; d = 1,op_nstepsDensityRO (indexRO)       (in # ions/cm^3)
!!
!!           op_tableTemperaturePA (t,indexPA) ; t = 1,op_nstepsTemperaturePA (indexPA)   (in eV, 1eV = 11405K)
!!           op_tableTemperaturePE (t,indexPE) ; t = 1,op_nstepsTemperaturePA (indexPE)   (in eV, 1eV = 11405K)
!!           op_tableTemperatureRO (t,indexRO) ; t = 1,op_nstepsTemperaturePA (indexRO)   (in eV, 1eV = 11405K)
!!
!!  The energy group boundaries will be stored in the 1-dimensional array:
!!
!!           op_energyGroupBoundaries (g) ; g = 1,op_ngroupEnergy+1     (in eV)
!!
!!  where
!!           op_energyGroupBoundaries (g)   = lower boundary of group 'g'
!!           op_energyGroupBoundaries (g+1) = upper boundary of group 'g'
!!
!! ARGUMENTS
!!
!!  tableName   : the name of the IONMIX file
!!  needPATable : if yes, Planck Absorption Opacities are needed from the IONMIX table
!!  needPETable : if yes, Planck   Emission Opacities are needed from the IONMIX table
!!  needROTable : if yes,         Rosseland Opacities are needed from the IONMIX table
!!  indexPA     : table counting index where Planck Absorption Opacities will be placed
!!  indexPE     : table counting index where Planck   Emission Opacities will be placed
!!  indexRO     : table counting index where Planck  Transport Opacities will be placed
!!
!!***
subroutine op_readIonmixTables (tableName,   &
                                needPATable, &
                                needPETable, &
                                needROTable, &
                                indexPA,     &
                                indexPE,     &
                                indexRO      )
  implicit none

  character (len=80), intent (in) :: tableName
  logical,            intent (in) :: needPATable
  logical,            intent (in) :: needPETable
  logical,            intent (in) :: needROTable
  integer,            intent (in) :: indexPA
  integer,            intent (in) :: indexPE
  integer,            intent (in) :: indexRO

  return
end subroutine op_readIonmixTables
