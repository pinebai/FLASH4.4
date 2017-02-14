!!****if* source/physics/Eos/EosMain/Tabulated/Hdf5TableRead/eos_tabReadOpacplotTables
!!
!! NAME
!!
!!  eos_tabReadOpacplotTables
!!
!! SYNOPSIS
!!
!!  call eos_tabReadOpacplotTables (character (in) :: tableName (len=80),
!!                            character   (in) :: groupName (len=80),
!!                            logical   (in) :: needZFTable,
!!                            logical   (in) :: needENTable,
!!                            logical   (in) :: needHCTable)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an OPACPLOT datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays: ! DEV: removed !
!!
!!             eos_tabIonizationTables  (t,d,g,indexZF)  (in cm^2/g)
!!             eos_tabIntEnergyTables    (t,d,g,indexEN)  (in cm^2/g)
!!             eos_tabHeatCpTables        (t,d,g,indexHC)  (in cm^2/g)
!!
!!    where:   t = temperature (eV) index
!!             d = ion number density (# ions/cm^3) index
!!             g = energy group (eV) index
!!       indexXX = table counting index for Opacity kind XX (XX = ZF,EN,HC)
!!
!!  The number of temperature and density indices will be stored in the 2-dimensional arrays: 
!!
!!         op_nstepsDensity (EOS_TABULAR_Z,indexZF)  !!DEV: outdated! - KW
!!         op_nstepsDensity (EOS_TABULAR_E,indexEN)
!!         op_nstepsDensity (EOS_TABULAR_C,indexHC)
!!
!!         op_nstepsTemperature (EOS_TABULAR_Z,indexZF)
!!         op_nstepsTemperature (EOS_TABULAR_E,indexEN)
!!         op_nstepsTemperature (EOS_TABULAR_C,indexHC)
!!
!!  The actual temperatures and densities will be stored in the 2-dimensional arrays: !!DEV: removed !
!!
!!         op_plasmaDensityZF     (d,indexZF) ; d = 1,op_nstepsDensity (EOS_TABULAR_Z,indexZF)       (in # ions/cm^3)
!!         op_plasmaDensityEN     (d,indexEN) ; d = 1,op_nstepsDensity (EOS_TABULAR_E,indexEN)       (in # ions/cm^3)
!!         op_plasmaDensityHC     (d,indexHC) ; d = 1,op_nstepsDensity (EOS_TABULAR_C,indexHC)       (in # ions/cm^3)
!!
!!         op_plasmaTemperatureZF (t,indexZF) ; t = 1,op_nstepsTemperature (EOS_TABULAR_Z,indexZF)   (in eV, 1eV = 11405K)
!!         op_plasmaTemperatureEN (t,indexEN) ; t = 1,op_nstepsTemperature (EOS_TABULAR_E,indexEN)   (in eV, 1eV = 11405K)
!!         op_plasmaTemperatureHC (t,indexHC) ; t = 1,op_nstepsTemperature (EOS_TABULAR_C,indexHC)   (in eV, 1eV = 11405K)
!!
!!
!! ARGUMENTS
!!
!!  tableName   : the name of the OPACPLOT file
!!  groupName   : the name of the OPACPLOT group (typically material name)
!!  indexZF     : table counting index where average ionization data will be placed
!!  indexEN     : table counting index where internal energy data will be placed
!!  indexHC     : table counting index where specific heat data will be placed
!!
!!***

#include "constants.h"

subroutine eos_tabReadOpacplotTables (tableName,   &
                                groupName, &
                                wanted, &
                                td, &
                                tbZF,tbEN,tbPR,tbHC,tbEntr)

  use Eos_data,         ONLY :  eos_smallT
  use eos_tabData,      ONLY :  eos_useLogTables,      &
                                EOS_TAB_NCOMP,         &
                                EOS_TAB_NALLTAB,       &
                                EOS_TAB_FOR_ION,       &
                                EOS_TAB_FOR_ELE,       &
                                EOS_TAB_FOR_MAT,       &
                                EOS_TABVT_ZF,          &
                                EOS_TABVT_EN,          &
                                EOS_TABVT_PR,          &
                                EOS_TABVT_HC,          &
                                EOS_TABVT_ENTR,        &
                                eosT_tableGroupDescT,  &
                                eosT_oneVarTablePT

  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  type(eosT_tableGroupDescT),intent(inout) :: td(:)
  type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr
  character (len=80), intent (in) :: tableName
  character (len=80), intent (in) :: groupName

  character (len=80) :: dummyLine
  character (len=81) :: nullTermTableName, nullTermGroupName


!!  needZFTable : if yes, average ionization data are needed from the OPACPLOT table
!!  needENTable : if yes, internal energy data are needed from the OPACPLOT table
!!  needHCTable : if yes,        specific heat data are needed from the OPACPLOT table
  logical :: needZFTable
  logical :: needENTables
  logical :: needPRTables
  logical :: needHCTables
  logical :: needEntrTables
  logical :: fileExists, doread

  integer :: ntemp, ndens 
  real, allocatable :: temperatures(:)
  real, allocatable :: densities(:)
  integer :: fileUnit
  integer :: notneededData
  integer :: step
  integer :: i,n,t,d,g
  integer :: ngroupsEnergy
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit

  real    :: dummyData
!!$  real    :: log10Density
  real    :: log10DensityStep
  real    :: log10DensityFirst
!!$  real    :: log10Temperature
  real    :: log10TemperatureStep
  real    :: log10TemperatureFirst
  real    :: maxTempStepFactor,maxDensStepFactor
  real    :: minTempStepFactor,minDensStepFactor

  ! Conversion from eV to K:
  real, parameter :: K = 11604.5221
  real, parameter :: joule2erg = 1.0e7

!   ...Check for existence
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[eos_tabReadOpacplotTables] ERROR: no Opacplot file found')
  end if

  nullTermTableName = trim(tableName)//char(0)
  nullTermGroupName = trim(groupName)//char(0)
  call eos_tabBrowsehdf5(nullTermTableName, nullTermGroupName, nstepsTemperature, nstepsDensity, ngroupsEnergy)

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[eos_tabReadOpacplotTables] ERROR: no OPACPLOT temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabReadOpacplotTables] ERROR: no OPACPLOT density grid found')
  end if

  if (ngroupsEnergy <= 0) then
      call Logfile_stampMessage('[eos_tabReadOpacplotTables] WARNING: no OPACPLOT energy group grid found')
  end if

  ntemp = nstepsTemperature
  ndens = nstepsDensity

  allocate(temperatures(ntemp))
  allocate(densities(ndens))
  allocate(tbZF(EOS_TAB_FOR_ELE)%table(ntemp,ndens)) ! avg ionization
  allocate(tbPR(EOS_TAB_FOR_ION)%table(ntemp,ndens)) ! ion pressure
  allocate(tbPR(EOS_TAB_FOR_ELE)%table(ntemp,ndens)) ! electron pressure
  allocate(tbEN(EOS_TAB_FOR_ION)%table(ntemp,ndens)) ! ion internal energy
  allocate(tbEN(EOS_TAB_FOR_ELE)%table(ntemp,ndens)) ! electron internal energy
  allocate(tbHC(EOS_TAB_FOR_ION)%table(ntemp,ndens)) ! ion heat capacity
  allocate(tbHC(EOS_TAB_FOR_ELE)%table(ntemp,ndens)) ! electron heat capacity
  allocate(tbEntr(EOS_TAB_FOR_ELE)%table(ntemp,ndens)) ! entropy

  call eos_read_hdf5(nullTermTableName, &
                     nullTermGroupName, &
                     temperatures, &
                     densities, &
                     tbZF(EOS_TAB_FOR_ELE)%table(ntemp,ndens), & !! avg ionization
                     tbPR(EOS_TAB_FOR_ION)%table(ntemp,ndens), & !! ion pressure
                     tbPR(EOS_TAB_FOR_ELE)%table(ntemp,ndens), & !! electron pressure
                     tbEN(EOS_TAB_FOR_ION)%table(ntemp,ndens), & !! ion internal energy
                     tbEN(EOS_TAB_FOR_ELE)%table(ntemp,ndens), & !! electron internal energy
                     tbHC(EOS_TAB_FOR_ION)%table(ntemp,ndens), & !! ion heat capacity
                     tbHC(EOS_TAB_FOR_ELE)%table(ntemp,ndens), & !! electron heat capacity
                     tbEntr(EOS_TAB_FOR_ELE)%table(ntemp,ndens))  !! entropy

  ! Convert the temperatures into K from eV:
  temperatures(:) = temperatures(:) * K
  log10TemperatureFirst = log10(Temperatures(1))
  minTempStepFactor = HUGE(minTempStepFactor)
  maxTempStepFactor = 1
  do t = 1,ntemp-1
     if(temperatures(t) < eos_smallT) then
        if(temperatures(t+1) > eos_smallT) &
             minTempStepFactor = min(minTempStepFactor,temperatures(t+1)/eos_smallT)
        maxTempStepFactor = max(maxTempStepFactor,temperatures(t+1)/eos_smallT)
     else
        minTempStepFactor = min(minTempStepFactor,temperatures(t+1)/temperatures(t))
        maxTempStepFactor = max(maxTempStepFactor,temperatures(t+1)/temperatures(t))
     end if
  end do

  log10DensityFirst = log10(densities(1))
  minDensStepFactor = HUGE(minDensStepFactor)
  maxDensStepFactor = 1
  do d = 1,ndens-1
     minDensStepFactor = min(maxDensStepFactor,densities(d+1)/densities(d))
     maxDensStepFactor = max(maxDensStepFactor,densities(d+1)/densities(d))
  end do

  needZFTable = wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF)
  needENTables = ANY(wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_EN))
  needPRTables = ANY(wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_PR))
  needHCTables = ANY(wanted(:,EOS_TABVT_HC))
  needEntrTables = ANY(wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_ENTR))

  if (needZFTable) then
      td(EOS_TABVT_ZF)%ntemp = nstepsTemperature
      td(EOS_TABVT_ZF)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_ZF)%Temperatures)) then
         allocate(td(EOS_TABVT_ZF)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_ZF)%Densities)) then
         allocate(td(EOS_TABVT_ZF)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_ZF)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_ZF)%Temperatures(step) = temperatures(step)
      end do

!!$     print*,'LBOUND(tbZF,1):',LBOUND(tbZF,1)
!!$     print*,'UBOUND(tbZF,1):',UBOUND(tbZF,1)

      if (td(EOS_TABVT_ZF)%isLog) then
         td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) )
         tbZF(EOS_TAB_FOR_ELE)%isLogData = .TRUE.
      end if
      if (tbZF(EOS_TAB_FOR_ELE)%isLogData) then
         tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity) = &
              log10(tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity))
      end if
  end if


!!!!!!!!!!!!!! Pressure tables !!!!!!!!!!!!!!!!!!!
  if (needPRTables) then
      td(EOS_TABVT_PR)%ntemp = nstepsTemperature
      td(EOS_TABVT_PR)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_PR)%Temperatures)) then
         allocate(td(EOS_TABVT_PR)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_PR)%Densities)) then
         allocate(td(EOS_TABVT_PR)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_PR)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_PR)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_PR)%isLog) then
         td(EOS_TABVT_PR)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_PR)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_PR)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_PR)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_PR)) tbPR(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_PR)
            if (doread) then
               tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                        tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg
               if (tbPR(n)%isLogData) then
                  tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                      log10(tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity))
               end if
            end if
      end do
  end if



!!!!!!!!!!!!!!! Energy Tables !!!!!!!!!

  if (needENTables) then
      td(EOS_TABVT_EN)%ntemp = nstepsTemperature
      td(EOS_TABVT_EN)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_EN)%Temperatures)) then
         allocate(td(EOS_TABVT_EN)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_EN)%Densities)) then
         allocate(td(EOS_TABVT_EN)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_EN)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_EN)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_EN)%isLog) then
         td(EOS_TABVT_EN)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_EN)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_EN)) tbEN(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_EN)
         if (doread) then
            tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg
            if (tbEN(n)%isLogData) then
               tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         end if
      end do
  end if


!!!!!!!!!! Heat Capacity !!!!!!!!!!!!

  if (needHCTables) then

      td(EOS_TABVT_HC)%ntemp = nstepsTemperature
      td(EOS_TABVT_HC)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_HC)%Temperatures)) then
         allocate(td(EOS_TABVT_HC)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_HC)%Densities)) then
         allocate(td(EOS_TABVT_HC)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_HC)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_HC)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_HC)%isLog) then
         td(EOS_TABVT_HC)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_HC)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_HC)) tbHC(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_HC)
         if (doread) then
            tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg / K
            if (tbHC(n)%isLogData) then
               tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         end if
      end do
  end if




!!!!!!!!! Entropy !!!!!!!!!!!
  if (needEntrTables) then
      td(EOS_TABVT_ENTR)%ntemp = nstepsTemperature
      td(EOS_TABVT_ENTR)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_ENTR)%Temperatures)) then
         allocate(td(EOS_TABVT_ENTR)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_ENTR)%Densities)) then
         allocate(td(EOS_TABVT_ENTR)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_ENTR)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_ENTR)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_ENTR)%isLog) then
         td(EOS_TABVT_ENTR)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_ENTR)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_ENTR)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_ENTR)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_ENTR)) tbEntr(n)%isLogData = .TRUE.
         end do
      end if

      n = EOS_TAB_FOR_ELE
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_ENTR)
         if (doread) then
            tbEntr(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbEntr(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg / K
            if (tbEntr(n)%isLogData) then
               tbEntr(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbEntr(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         end if
  end if

  deallocate(temperatures)
  deallocate(densities)

  return

end subroutine eos_tabReadOpacplotTables
