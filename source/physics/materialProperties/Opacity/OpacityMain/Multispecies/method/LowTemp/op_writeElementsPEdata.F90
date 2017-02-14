!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_writeElementsPEdata
!!
!! NAME
!!
!!  op_writeElementsPEdata
!!
!! SYNOPSIS
!!
!!  call op_writeElementsPEdata ()
!!
!! DESCRIPTION
!!
!!  This routine writes all the current atomic elements data for determining the photoelectric cross
!!  sections to a specific file. The routine is meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeElementsPEdata ()

  use Opacity_data,     ONLY : op_atomName,             &
                               op_totalElements,        &
                               op_element2AtomicNumber

  use op_lowTempData,   ONLY : op_elementAij4,          &
                               op_elementJmax,          &
                               op_elementPEenergyRange

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

# include "Opacity.h"

  logical :: arrayExists
  logical :: goodShape

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: j,Z
  integer :: jmax
  integer :: n
!
!
!   ...Do the printout.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "Biggs_and_Lighthill_Elements_PE_data.txt", &
        form = 'formatted')

  write (fileUnit,*)
  write (fileUnit,*) '   Biggs & Lighthill Elements Photoelectron Cross Section Data'
  write (fileUnit,*) '   ==========================================================='
  write (fileUnit,*)

  do n = 1,op_totalElements

     Z = op_element2AtomicNumber (n)

     write (fileUnit,*)
     write (fileUnit,'(2X,I3,2X,A12)')       Z, op_atomName (Z)
     write (fileUnit,'(4X,A1,      6X,A5,3X,   3X,A6,2X,  4(1X,A9,1X))') &
                          'J',       'START',   'FINISH','A (I,J,1)','A (I,J,2)','A (I,J,3)','A (I,J,4)'
     write (fileUnit,*) '   ---------------------------------------------------------------------'

     jmax = op_elementJmax (n)

     do j = 1,jmax
        write (fileUnit,'(3X,I2,                          2(1X,F10.4,1X),            4(1X,ES10.3))') &
                              j,  op_elementPEenergyRange (LOW:HIGH,j,n),   op_elementAij4 (1:4,j,n)
     end do

  end do

  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeElementsPEdata
