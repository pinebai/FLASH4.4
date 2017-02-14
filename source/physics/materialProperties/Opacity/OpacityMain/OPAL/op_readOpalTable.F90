!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_readOpalTable
!!
!! NAME
!!
!!  op_readOpalTable
!!
!! SYNOPSIS
!!
!!  call op_readOpalTable (character (in) :: tableName (len=80),
!!                             logical   (in) :: needLowTable,
!!                             logical   (in) :: needHighTable,
!!                             logical   (in) :: needROTable)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an OPAL datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays:
!!
!!             op_PlanckAbsorptionTables (t,d,g,indexLOWT)  (in cm^2/g)
!!             op_PlanckEmissionTables   (t,d,g,indexHIGHT)  (in cm^2/g)
!!             op_RosselandTables        (t,d,g,indexRO)  (in cm^2/g)
!!
!!    where:   t = temperature (K) index
!!             d = ion number density (# ions/cm^3) index
!!             g = energy group (eV) index
!!       indexXX = table counting index for Opacity kind XX (XX = LowT,HighT,RO)
!!
!!  The number of temperature and density indices will be stored in the 1-dimensional arrays: 
!!
!!             op_nstepsDensityLowT  (indexLOWT)
!!             op_nstepsDensityHighT  (indexHIGHT)
!!             op_nstepsDensityRO  (indexRO)
!!
!!             op_nstepsTemperatureLowT  (indexLOWT)
!!             op_nstepsTemperatureHighT  (indexHIGHT)
!!             op_nstepsTemperatureRO  (indexRO)
!!
!!  The actual temperatures and densities will be stored in the 2-dimensional arrays:
!!
!!           op_tableDensityLowT     (d,indexLOWT) ; d = 1,op_nstepsDensityLowT (indexLOWT)       (in # ions/cm^3)
!!           op_tableDensityHighT     (d,indexHIGHT) ; d = 1,op_nstepsDensityHighT (indexHIGHT)       (in # ions/cm^3)
!!           op_tableDensityRO     (d,indexRO) ; d = 1,op_nstepsDensityRO (indexRO)       (in # ions/cm^3)
!!
!!           op_tableTemperatureLowT (t,indexLOWT) ; t = 1,op_nstepsTemperatureLowT (indexLOWT)   (in K)
!!           op_tableTemperatureHighT (t,indexHIGHT) ; t = 1,op_nstepsTemperatureLowT (indexHIGHT)   (in K)
!!           op_tableTemperatureRO (t,indexRO) ; t = 1,op_nstepsTemperatureLowT (indexRO)   (in K)
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
!!  tableName   : the name of the OPAL file
!!
!!
!!***
subroutine op_readOpalTable (tableName,   &
                                      td, &
                                      tb  )

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups
  use Opacity_data,  ONLY : op_globalMe

  use op_opalData, ONLY :  opT_tableGroupDescT,  &
                                opT_oneVarTablePT
  use op_opalData,  ONLY : op_useLogTables

  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Opacity.h"

  character (len=80), intent (in) :: tableName
       type(opT_tableGroupDescT),intent(inout) :: td
       type(opT_oneVarTablePT),pointer :: tb

  real, parameter :: ten   = 10.0

  character (len=20) :: fmtStr
  character (len=120) :: line1
  integer :: form, version, nstepsR
  real    :: X,Z, logRmin,logRmax, logTmin,logTmax

  logical :: fileExists

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  integer :: t,d
  integer,parameter :: nEnergyGroups = 1
  integer :: nstepsDensity
  integer :: nstepsTemperature

  integer :: ntemp, ndens
  real, allocatable :: logRs(:)
  real, allocatable :: logTs(:)
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[op_readOpalTable] ERROR: no OPAL file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature, density and energy group grids. Abort the calculation,
!      if any of the grids is not found.
!
!

!
!
!   ...Read the temperature and density grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  if (op_globalMe == MASTER_PE) print*,'R OPAL file ', trim(tableName),':'
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,'R',line1
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,'R',line1
!    form     version        X           Z          logRs    logR_min    logR_max       logTs    logT_min    logT_max
!       1          37    0.000000    0.020000          46   -8.000000    1.000000         105    2.700000    4.500000
!88format(I8,1x,I11,2G)
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,'R',line1
  read (line1,*) form,version,X,Z, nstepsR, logRmin,logRmax, nstepsTemperature,logTmin,logTmax
!!$  if (op_globalMe == MASTER_PE) print*, form,version,X,Z, nstepsR, logRmin,logRmax, nstepsTemperature,logTmin,logTmax
  nstepsDensity = nstepsR
  ntemp = nstepsTemperature
  ndens = nstepsDensity
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,'R',line1
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,'R',line1

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_readOpalTable] ERROR: no OPAL temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_readOpalTable] ERROR: no OPAL density grid found')
  end if

  allocate(logRs(nstepsDensity))
  allocate(logTs(nstepsTemperature))

70 format('(8x,',i4,'(e8.3))')
   write(fmtStr,70) nstepsDensity
#ifdef DEBUG_OPACITY
   print*,'fmtStr:',fmtStr
#endif

  read (fileUnit,fmtStr) (logRs(d),d=1,nstepsDensity)
  read (fileUnit,'(A120)') line1

!
!
!   ...Establish the temperature and density grid for the Rosseland case (if needed)
!      and read in the opacities.
!
!

  if (.TRUE.) then

      td%ntemp = nstepsTemperature
      td%ndens = nstepsDensity
      if (.NOT.associated(td%Temperatures)) then
         allocate(td%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td%Densities)) then
         allocate(td%Densities(ndens))
      end if

      allocate(tb%table(nstepsTemperature,nstepsDensity))
#ifdef DEBUG_OPACITY
      print*,'LBOUND(tb%table,1):',LBOUND(tb%table,1)
      print*,'UBOUND(tb%table,1):',UBOUND(tb%table,1)
#endif

80    format('((',i4,'F8.3))')
      write(fmtStr,80) nstepsDensity+1
#ifdef DEBUG_OPACITY
      print*,'fmtStr:',fmtStr
#endif
!!$      read (fileUnit,fmtStr) (logTs(t), (tb%table(t,d), d = 1,nstepsDensity), t = 1,nstepsTemperature    )
      do t=1,nstepsTemperature
         read (fileUnit,fmtStr) logTs(t), (tb%table(t,d), d = 1,nstepsDensity)
#ifdef DEBUG_OPACITY
         write(*       ,fmtStr) logTs(t), (tb%table(t,d), d = 1,nstepsDensity)
#endif
      end do


      if (td%isLog) then
         td%Densities    (1:nstepsDensity)     =       logRs  (1:nstepsDensity)
         td%Temperatures (1:nstepsTemperature) =       logTs  (1:nstepsTemperature)
         tb%isLogData = .TRUE.
      else
         td%Densities    (1:nstepsDensity)     = ten**(logRs  (1:nstepsDensity) )
         td%Temperatures (1:nstepsTemperature) = ten**(logTs  (1:nstepsTemperature) )
      end if
      if (.NOT. tb%isLogData) then
         tb%table(1:nstepsTemperature,1:nstepsDensity) = &
              ten**(tb%table(1:nstepsTemperature,1:nstepsDensity))
      end if
  end if
!
!
!   ...Close the OPAL file and deallocate temperature/density arrays
!
!
  read (fileUnit,'(A120)') line1
!#ifdef DEBUG_OPACITY
  if (op_globalMe == MASTER_PE) print*,'R',trim(line1),'.'
!#endif
  deallocate(logRs)
  deallocate(logTs)
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readOpalTable
