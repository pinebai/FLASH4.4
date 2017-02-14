!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_writeCompletePEdata
!!
!! NAME
!!
!!  op_writeCompletePEdata
!!
!! SYNOPSIS
!!
!!  call op_writeCompletePEdata ()
!!
!! DESCRIPTION
!!
!!  This routine writes all the data for determining the photoelectric cross sections to a specific file.
!!  The routine is meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeCompletePEdata ()

  use Opacity_data,     ONLY : op_atomName

  use op_lowTempData,   ONLY : op_maxElements,      &
                               op_Aij4,             &
                               op_Jmax,             &
                               op_PEenergyRange

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

# include "Opacity.h"

  logical :: arrayExists
  logical :: goodShape

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: j,Z
  integer :: jmax
!
!
!   ...Do the printout.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "Biggs_and_Lighthill_PE_data.txt", &
        form = 'formatted')

  write (fileUnit,*)
  write (fileUnit,*) '   Biggs & Lighthill Photoelectron Cross Section Data'
  write (fileUnit,*) '   =================================================='
  write (fileUnit,*)

  do Z = 1,op_maxElements

     write (fileUnit,*)
     write (fileUnit,'(2X,I3,2X,A12)')       Z, op_atomName (Z)
     write (fileUnit,'(4X,A1,      6X,A5,3X,   3X,A6,2X,  4(1X,A9,1X))') &
                          'J',       'START',   'FINISH','A (I,J,1)','A (I,J,2)','A (I,J,3)','A (I,J,4)'
     write (fileUnit,*) '   ---------------------------------------------------------------------'

     jmax = op_Jmax (Z)

     do j = 1,jmax
        write (fileUnit,'(3X,I2,                   2(1X,F10.4,1X),        4(1X,ES10.3))') &
                              j,  op_PEenergyRange (LOW:HIGH,j,Z),   op_Aij4 (1:4,j,Z)
     end do

  end do

  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeCompletePEdata
