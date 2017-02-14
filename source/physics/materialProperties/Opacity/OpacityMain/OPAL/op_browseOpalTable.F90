!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_browseOpalTable
!!
!! NAME
!!
!!  op_browseOpalTable
!!
!! SYNOPSIS
!!
!!  call op_browseOpalTable (character (in)  :: tableName (len=80),
!!                              integer   (out) :: nstepsDensity,
!!                              integer   (out) :: nstepsTemperature)
!!
!! DESCRIPTION
!!
!!  This routine browses through the tabulated opacities from an OPAL datafile output in
!!  order to extract the number of steps for both the density and the temperature grid with
!!  which the OPAL tables were generated.
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the OPAL file
!!  nstepsDensity     : the size of the         Rosseland density grid returned
!!  nstepsTemperature : the size of the         Rosseland temperature grid returned
!!
!!***

#include "constants.h"

subroutine op_browseOpalTable (tableName,                        &
                                               nstepsDensity,     &
                                               nstepsTemperature  )

  use Driver_interface,  ONLY : Driver_abortFlash
  use Opacity_data,  ONLY : op_globalMe

  implicit none

  character (len=80), intent (in)  :: tableName
  integer,            intent (out) :: nstepsDensity
  integer,            intent (out) :: nstepsTemperature

  character (len=120) :: line1
  integer :: form, version, nstepsR
  real    :: X,Z, logRmin,logRmax, logTmin,logTmax
  logical :: fileExists

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
!
!
!   ...Check and open the OPAL opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
     if (op_globalMe==MASTER_PE) &
          print*,'[op_browseOpalTable] ERROR: OPAL file not found: ',tableName 
       call Driver_abortFlash ('[op_browseOpalTable] ERROR: no OPAL file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature and density grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  if (op_globalMe == MASTER_PE) print*,'OPAL file ', trim(tableName),':'
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,line1
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,line1
!    form     version        X           Z          logRs    logR_min    logR_max       logTs    logT_min    logT_max
!       1          37    0.000000    0.020000          46   -8.000000    1.000000         105    2.700000    4.500000
!88format(I8,1x,I11,2G)
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,line1
  read (line1,*) form,version,X,Z, nstepsR, logRmin,logRmax, nstepsTemperature,logTmin,logTmax
!!$  if (op_globalMe == MASTER_PE) print*, form,version,X,Z, nstepsR, logRmin,logRmax, nstepsTemperature,logTmin,logTmax
  nstepsDensity = nstepsR
  read (fileUnit,'(A120)') line1
  if (op_globalMe == MASTER_PE) print*,line1
  read (fileUnit,'(A120)') line1
!!$  if (op_globalMe == MASTER_PE) print*

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_browseOpalTable] ERROR: no OPAL temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_browseOpalTable] ERROR: no OPAL density grid found')
  end if
!
!
!
!   ...Close the OPAL file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_browseOpalTable
