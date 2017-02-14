!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_readIonmixTables
!!
!! NAME
!!
!!  op_readIonmixTables
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
!!    where:   t = temperature (K) index
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
!!           op_tableTemperaturePA (t,indexPA) ; t = 1,op_nstepsTemperaturePA (indexPA)   (in K)
!!           op_tableTemperaturePE (t,indexPE) ; t = 1,op_nstepsTemperaturePA (indexPE)   (in K)
!!           op_tableTemperatureRO (t,indexRO) ; t = 1,op_nstepsTemperaturePA (indexRO)   (in K)
!!
!!  The energy group boundaries will be stored temporarily in the 1-dimensional array:
!!
!!           op_tabulatedEnergyBoundaries (g) ; g = 1,op_nEnergyGroups+1     (in eV)
!!
!!  where
!!           op_tabulatedEnergyBoundaries (g)   = lower boundary of group 'g'
!!           op_tabulatedEnergyBoundaries (g+1) = upper boundary of group 'g'
!!
!!  and will be checked against the energy group boundaries stored for the opacity unit.
!!  A discrepancy within the predefined energy difference tolerance will signal an inconsistency
!!  in the tables generated for the current problem.
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
!! NOTES
!!
!!  Since the IONMIX tables are produced with temperature units in eV, we must convert these
!!  here into units of K. If the use of logarithmic Tables has been specified, we need to
!!  be careful to substitute the exact zeros of the tables by the smallest representable
!!  positive number in order to avoid NaN's when taking the logarithms.
!!
!!***
subroutine op_readIonmixTables (tableName,   &
                                needPATable, &
                                needPETable, &
                                needROTable, &
                                indexPA,     &
                                indexPE,     &
                                indexRO      )

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups,             &
                                op_energyGroupBoundaries

  use op_tabulatedData,  ONLY : op_useLogTables,              &
                                op_tableEnergyTolerance,      &
                                op_nstepsDensityPA,           &
                                op_nstepsDensityPE,           &
                                op_nstepsDensityRO,           &
                                op_nstepsTemperaturePA,       &
                                op_nstepsTemperaturePE,       &
                                op_nstepsTemperatureRO,       &
                                op_tabulatedEnergyBoundaries, &
                                op_tableDensityPA,            &
                                op_tableDensityPE,            &
                                op_tableDensityRO,            &
                                op_tableTemperaturePA,        &
                                op_tableTemperaturePE,        &
                                op_tableTemperatureRO,        &
                                op_RosselandTables,           &
                                op_PlanckAbsorptionTables,    &
                                op_PlanckEmissionTables

  use op_numericsData,   ONLY : zero,ten,                     &
                                op_eV2Kelvin,                 &
                                op_smallestPositiveNumber

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Opacity.h"

  character (len=80), intent (in) :: tableName
  logical,            intent (in) :: needPATable
  logical,            intent (in) :: needPETable
  logical,            intent (in) :: needROTable
  integer,            intent (in) :: indexPA
  integer,            intent (in) :: indexPE
  integer,            intent (in) :: indexRO

  character (len=80) :: dummyLine

  logical :: fileExists

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  integer :: notneededData
  integer :: step
  integer :: n,t,d,g
  integer :: nEnergyGroups
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: skipTablePA
  integer :: skipTableRO

  real    :: dummyData
  real    :: energyDifference
  real    :: log10Density
  real    :: log10DensityStep
  real    :: log10DensityFirst
  real    :: log10Temperature
  real    :: log10TemperatureStep
  real    :: log10TemperatureFirst
  real    :: newZero
  real    :: temperatureFirst
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[op_readIonmixTables] ERROR: no IONMIX file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature, density and energy group grids. Abort the calculation,
!      if any of the grids is not found. Check also, if the number of energy groups
!      read from the file corresponds to the one of the current run.
!
!
  read (fileUnit,'(2I10)') nstepsTemperature , nstepsDensity
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine

  read (fileUnit,'(4E12.6,I12)') log10DensityStep,       &
                                 log10DensityFirst,      &
                                 log10TemperatureStep,   &
                                 log10TemperatureFirst,  &
                                 nEnergyGroups

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_readIonmixTables] ERROR: no IONMIX temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_readIonmixTables] ERROR: no IONMIX density grid found')
  end if

  if (nEnergyGroups <= 0) then
      call Driver_abortFlash ('[op_readIonmixTables] ERROR: no IONMIX energy group grid found')
  end if

  if (nEnergyGroups /= op_nEnergyGroups) then
      call Driver_abortFlash ('[op_readIonmixTables] ERROR: bad size of IONMIX energy group grid')
  end if
!
!
!   ...Perform the temperature eV -> K conversion. This needs to be done only on the
!      first temperature of the grid. The temperature step only defines the multiplicative
!      factor as we move along the temperature grid.
!
!
  temperatureFirst      = ten ** log10TemperatureFirst             !  log(T) -> T in eV
  temperatureFirst      = temperatureFirst * op_eV2Kelvin          ! T in eV -> T in K
  log10TemperatureFirst = log10 (temperatureFirst)                 !       T -> log(T)
!
!
!   ...Skip not needed data from the IONMIX file.
!
!
  notneededData =   nstepsDensity * nstepsTemperature

  read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
!
!
!   ...Read the energy group boundaries from the IONMIX file and compare with those stored.
!
!
  read (fileUnit,'(4E12.6)') (op_tabulatedEnergyBoundaries (g), g = 1,nEnergyGroups + 1)

  do g = 1,nEnergyGroups+1

     energyDifference = abs (op_tabulatedEnergyBoundaries (g) - op_energyGroupBoundaries (g)) / &
          op_energyGroupBoundaries(g)

     if (energyDifference > op_tableEnergyTolerance) then
         call Driver_abortFlash ('[op_readIonmixTables] ERROR: IONMIX / Opacity energy group mismatch')
     end if
  end do
!
!
!   ...Establish the temperature and density grid for the Rosseland case (if needed)
!      and read in the opacities.
!
!
  skipTableRO = 1
  skipTablePA = 1

  if (needROTable) then

      op_nstepsDensityRO     (indexRO) = nstepsDensity
      op_nstepsTemperatureRO (indexRO) = nstepsTemperature

      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         op_tableDensityRO (step,indexRO) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         op_tableTemperatureRO (step,indexRO) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

      read (fileUnit,'(4E12.6)') (((op_RosselandTables (t,d,g,indexRO), t = 1,nstepsTemperature) &
                                                                      , d = 1,nstepsDensity    ) &
                                                                      , g = 1,nEnergyGroups    )
      if (op_useLogTables) then
                   op_tableDensityRO      (1:nstepsDensity,    indexRO) &
          = log10 (op_tableDensityRO      (1:nstepsDensity,    indexRO) )
                   op_tableTemperatureRO  (1:nstepsTemperature,indexRO) &
          = log10 (op_tableTemperatureRO  (1:nstepsTemperature,indexRO) )

          newZero = op_smallestPositiveNumber

                   op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) &
          =   max (op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO),newZero)
                   op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) &
          = log10 (op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) )
      end if

      skipTableRO = 0

  end if
!
!
!   ...Establish the temperature and density grid for the Planck Absorption case (if needed)
!      and read in the opacities.
!
!
  if (needPATable) then

      op_nstepsDensityPA     (indexPA) = nstepsDensity
      op_nstepsTemperaturePA (indexPA) = nstepsTemperature

      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         op_tableDensityPA (step,indexPA) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         op_tableTemperaturePA (step,indexPA) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

      notneededData = skipTableRO * nstepsDensity * nstepsTemperature * nEnergyGroups

      if (notneededData > 0) then
          read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
      end if

      read (fileUnit,'(4E12.6)') (((op_PlanckAbsorptionTables (t,d,g,indexPA), t = 1,nstepsTemperature) &
                                                                             , d = 1,nstepsDensity    ) &
                                                                             , g = 1,nEnergyGroups    )
      if (op_useLogTables) then
                   op_tableDensityPA         (1:nstepsDensity,    indexPA) &
          = log10 (op_tableDensityPA         (1:nstepsDensity,    indexPA) )
                   op_tableTemperaturePA     (1:nstepsTemperature,indexPA) &
          = log10 (op_tableTemperaturePA     (1:nstepsTemperature,indexPA) )

          newZero = op_smallestPositiveNumber

                   op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) &
          =   max (op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA),newZero)
                   op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) &
          = log10 (op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) )
      end if

      skipTableRO = 0
      skipTablePA = 0

  end if
!
!
!   ...Establish the temperature and density grid for the Planck Emission case (if needed)
!      and read in the opacities.
!
!
  if (needPETable) then

      op_nstepsDensityPE     (indexPE) = nstepsDensity
      op_nstepsTemperaturePE (indexPE) = nstepsTemperature

      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         op_tableDensityPE (step,indexPE) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         op_tableTemperaturePE (step,indexPE) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

      notneededData =   skipTableRO * nstepsDensity * nstepsTemperature * nEnergyGroups &
                      + skipTablePA * nstepsDensity * nstepsTemperature * nEnergyGroups

      if (notneededData > 0) then
           read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
      end if

      read (fileUnit,'(4E12.6)') (((op_PlanckEmissionTables (t,d,g,indexPE), t = 1,nstepsTemperature) &
                                                                           , d = 1,nstepsDensity    ) &
                                                                           , g = 1,nEnergyGroups    )
      if (op_useLogTables) then
                   op_tableDensityPE       (1:nstepsDensity,    indexPE) &
          = log10 (op_tableDensityPE       (1:nstepsDensity,    indexPE) )
                   op_tableTemperaturePE   (1:nstepsTemperature,indexPE) &
          = log10 (op_tableTemperaturePE   (1:nstepsTemperature,indexPE) )

          newZero = op_smallestPositiveNumber

                   op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) &
          =   max (op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE),newZero)
                   op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) &
          = log10 (op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) )
      end if

  end if
!
!
!   ...Close the IONMIX file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readIonmixTables
