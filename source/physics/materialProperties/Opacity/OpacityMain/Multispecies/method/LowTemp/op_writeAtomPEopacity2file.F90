!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_writeAtomPEopacity2file
!!
!! NAME
!!
!!  op_writeAtomPEopacity2file
!!
!! SYNOPSIS
!!
!!  call op_writeAtomPEopacity2file (integer (in) :: Z,
!!                                   real    (in) :: Elower,
!!                                   real    (in) :: Eupper,
!!                                   integer (in) :: nPoints)
!!
!! DESCRIPTION
!!
!!  This routine writes the photoelectron opacities for element Z between the
!!  energy range Elower - Eupper in nPoints points to a specific file.
!!  The routine is meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!  Z       : atomic number
!!  Elower  : lower energy bound (in keV)
!!  Eupper  : upper energy bound (in keV)
!!  nPoints : number of points to be written to file
!!
!!***
subroutine op_writeAtomPEopacity2file (Z, Elower, Eupper, nPoints)

  use op_lowTempData,   ONLY : op_Aij4,           &
                               op_Jmax,           &
                               op_PEenergyRange

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

# include "Opacity.h"

  integer, intent (in) :: nPoints
  integer, intent (in) :: Z
  real,    intent (in) :: Elower
  real,    intent (in) :: Eupper

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: i,j,n
  integer :: jmax
  integer :: nSteps

  real    :: Elow,Ehigh
  real    :: Energy
  real    :: Estep
  real    :: Opacity
!
!
!   ...Check, if the PE opacity data is still available.
!
!
  if (     (.not.allocated (op_Aij4))          &
      .or. (.not.allocated (op_Jmax))          &
      .or. (.not.allocated (op_PEenergyRange)) ) then
       call Driver_abortFlash ('[op_writeAtomPEopacity2file] ERROR: No data available')
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "Atom_opacities_plot.txt", &
        form = 'formatted')
!
!
!   ...Write the atomic PE data.
!
!
  jmax = op_Jmax (Z)

  nSteps = nPoints - 1
  Estep  = (Eupper - Elower) / real (nSteps)

  do n = 1,nPoints

     Energy = Elower + real (n - 1) * Estep

     write (*,*) ' point #, energy = ',n,Energy


     do j = 1,jmax

        Elow  = op_PEenergyRange (LOW ,j,Z)
        Ehigh = op_PEenergyRange (HIGH,j,Z)

        write (*,*) ' Elow, Ehigh = ',Elow, Ehigh

        if ((Elow <= Energy) .and. (Energy < Ehigh)) then

             Opacity = 0.0
             do i = 1,4
                Opacity = Opacity + op_Aij4 (i,j,Z) / (Energy ** i)
             end do

             write (fileUnit,'(2ES10.3)') Energy,Opacity

             exit

        end if

     end do
  end do

  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeAtomPEopacity2file
