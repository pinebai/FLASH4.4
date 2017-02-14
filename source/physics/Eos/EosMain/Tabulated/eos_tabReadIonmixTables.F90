!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabReadIonmixTables
!!
!! NAME
!!
!!  eos_tabReadIonmixTables
!!
!! SYNOPSIS
!!
!!  call eos_tabReadIonmixTables (character (in) :: tableName (len=80),
!!                            logical   (in) :: needZFTable,
!!                            logical   (in) :: needENTable,
!!                            logical   (in) :: needHCTable)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an IONMIX datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays: !DEV: removed !
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
!!  tableName   : the name of the IONMIX file
!!  needZFTable : if yes, average ionization data are needed from the IONMIX table
!!  needENTable : if yes, internal energy data are needed from the IONMIX table
!!  needHCTable : if yes,        specific heat data are needed from the IONMIX table
!!  indexZF     : table counting index where average ionization data will be placed
!!  indexEN     : table counting index where internal energy data will be placed
!!  indexHC     : table counting index where specific heat data will be placed
!!
!!***
subroutine eos_tabReadIonmixTables (tableName,   &
                                wanted, &
                                td,&
                                tbZF,tbEN,tbHC)

  use eos_tabData,      ONLY : ten,                      &
                                eos_useLogTables,          &
                                EOS_TAB_NCOMP,         &
                                EOS_TAB_NALLTAB,         &
                                EOS_TAB_FOR_ION,         &
                                EOS_TAB_FOR_ELE,         &
                                EOS_TAB_FOR_MAT,         &
                                EOS_TABVT_ZF,         &
                                EOS_TABVT_EN,         &
                                EOS_TABVT_HC,         &
                                eosT_tableGroupDescT, &
                                eosT_oneVarTablePT

  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  type(eosT_tableGroupDescT),intent(inout) :: td(:)
  type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbHC
  character (len=80), intent (in) :: tableName

  character (len=80) :: dummyLine

  logical :: needZFTable
  logical :: needENTable
  logical :: needHCTable
  logical :: fileExists

  integer :: ntemp, ndens
  integer :: fileUnit
  integer :: notneededData
  integer :: step
  integer :: n,t,d,g
  integer :: ngroupsEnergy
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit

  real    :: dummyData
  real    :: log10Density
  real    :: log10DensityStep
  real    :: log10DensityFirst
  real    :: log10Temperature
  real    :: log10TemperatureStep
  real    :: log10TemperatureFirst

  ! Conversion from eV to K:
  real, parameter :: K = 11604.5221
  real, parameter :: joule2erg = 1.0e7
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[eos_tabReadIonmixTables] ERROR: no IONMIX file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature, density and energy group grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  read (fileUnit,'(2I10)') nstepsTemperature , nstepsDensity
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine

  read (fileUnit,'(4E12.6,I12)') log10DensityStep,       &
                                 log10DensityFirst,      &
                                 log10TemperatureStep,   &
                                 log10TemperatureFirst,  &
                                 ngroupsEnergy

  ! Convert the first temperature into K from eV:
  log10TemperatureFirst = log10(10.0**log10TemperatureFirst * K)

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[eos_tabReadIonmixTables] ERROR: no IONMIX temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabReadIonmixTables] ERROR: no IONMIX density grid found')
  end if

  if (ngroupsEnergy <= 0) then
      call Logfile_stampMessage('[eos_tabReadIonmixTables] WARNING: no IONMIX energy group grid found')
  end if

  ntemp = nstepsTemperature
  ndens = nstepsDensity


  needZFTable = wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF)
  needENTable = wanted(EOS_TAB_FOR_MAT,EOS_TABVT_EN)
  needHCTable = ANY(wanted(:,EOS_TABVT_HC))
!
!
!   ...Establish the temperature and density grid for the average ionization (if needed)
!      and read in the data.
!
!

  if (needZFTable) then

      td(EOS_TABVT_ZF)%ntemp = nstepsTemperature
      td(EOS_TABVT_ZF)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_ZF)%Temperatures)) then
         allocate(td(EOS_TABVT_ZF)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_ZF)%Densities)) then
         allocate(td(EOS_TABVT_ZF)%Densities(ndens))
      end if

      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         td(EOS_TABVT_ZF)%Densities(step) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         td(EOS_TABVT_ZF)%Temperatures(step) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

!!$     print*,'LBOUND(tbZF,1):',LBOUND(tbZF,1)
!!$     print*,'UBOUND(tbZF,1):',UBOUND(tbZF,1)
      allocate(tbZF(EOS_TAB_FOR_ELE)%table(nstepsTemperature,nstepsDensity))
      read (fileUnit,'(4E12.6)') ((tbZF(EOS_TAB_FOR_ELE)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

      if (td(EOS_TABVT_ZF)%isLog) then
         td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) )
         tbZF(EOS_TAB_FOR_ELE)%isLogData = .TRUE.
      end if
      if (tbZF(EOS_TAB_FOR_ELE)%isLogData) then
         tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity) = &
              log10(tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity))
      end if

  else
      !   ...Skip not needed data from the IONMIX file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  end if
!
!
!   ...Establish the temperature and density grid for the internal energy (if needed)
!      and read in the data.
!
!
  if (needENTable) then

      td(EOS_TABVT_EN)%ntemp = nstepsTemperature
      td(EOS_TABVT_EN)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_EN)%Temperatures)) then
         allocate(td(EOS_TABVT_EN)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_EN)%Densities)) then
         allocate(td(EOS_TABVT_EN)%Densities(ndens))
      end if

      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         td(EOS_TABVT_EN)%Densities(step) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         td(EOS_TABVT_EN)%Temperatures(step) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

      allocate(tbEN(EOS_TAB_FOR_MAT)%table(nstepsTemperature,nstepsDensity))
      read (fileUnit,'(4E12.6)') ((tbEN(EOS_TAB_FOR_MAT)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

      if (td(EOS_TABVT_EN)%isLog) then
         td(EOS_TABVT_EN)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_EN)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) )
         tbEN(EOS_TAB_FOR_MAT)%isLogData = .TRUE.
      end if
      tbEN(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) = &
           tbEN(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg
      if (tbEN(EOS_TAB_FOR_MAT)%isLogData) then
         tbEN(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) = &
              log10(tbEN(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity))
      end if

  else
      !   ...Skip not needed data from the IONMIX file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  end if
!
!
!   ...Establish the temperature and density grid for the specific heat capacity (if needed)
!      and read in the data.
!
!
  if (needHCTable) then

      td(EOS_TABVT_HC)%ntemp = nstepsTemperature
      td(EOS_TABVT_HC)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_HC)%Temperatures)) then
         allocate(td(EOS_TABVT_HC)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_HC)%Densities)) then
         allocate(td(EOS_TABVT_HC)%Densities(ndens))
      end if


      log10Density = log10DensityFirst
      do step = 1,nstepsDensity
         td(EOS_TABVT_HC)%Densities(step) = ten ** log10Density
         log10Density = log10Density + log10DensityStep
      end do

      log10Temperature = log10TemperatureFirst
      do step = 1,nstepsTemperature
         td(EOS_TABVT_HC)%Temperatures(step) = ten ** log10Temperature
         log10Temperature = log10Temperature + log10TemperatureStep
      end do

      allocate(tbHC(EOS_TAB_FOR_MAT)%table(nstepsTemperature,nstepsDensity))
      read (fileUnit,'(4E12.6)') ((tbHC(EOS_TAB_FOR_MAT)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

      if (td(EOS_TABVT_HC)%isLog) then
         td(EOS_TABVT_HC)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_HC)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) )
         tbHC(EOS_TAB_FOR_MAT)%isLogData = .TRUE.
      end if
      tbHC(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) = &
           tbHC(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg / K
      if (tbHC(EOS_TAB_FOR_MAT)%isLogData) then
         tbHC(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity) = &
              log10(tbHC(EOS_TAB_FOR_MAT)%table(1:nstepsTemperature,1:nstepsDensity))
      end if

      !
      !   ...Stop reading here, we do not care about the remainder of the file.
      !
      !


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
end subroutine eos_tabReadIonmixTables
