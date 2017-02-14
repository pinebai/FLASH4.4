!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_setEnergyGroupBoundaries
!!
!! NAME
!!
!!  op_setEnergyGroupBoundaries
!!
!! SYNOPSIS
!!
!!  call op_setEnergyGroupBoundaries ()
!!
!! DESCRIPTION
!!
!!  This routine sets the energy group boundaries from the data residing in the RadTrans
!!  unit. This is done in order to compare the data files for the tabulated opacities and
!!  to catch any inconsistencies.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setEnergyGroupBoundaries ()

  use Opacity_data,         ONLY : op_nEnergyGroups,         &
                                   op_energyGroupBoundaries

  use op_numericsData,      ONLY : zero,      &
                                   op_erg2eV

  use RadTrans_interface,   ONLY : RadTrans_mgdGetBound
  use Driver_interface,     ONLY : Driver_abortFlash

  implicit none

  integer :: bound
  integer :: nTotalBounds

  real    :: Eerg,EeV
!
!
!    ...Loop over all the boundaries and set the energies in eV.
!
!
  nTotalBounds = op_nEnergyGroups + 1

  do bound = 1,nTotalBounds

     call RadTrans_mgdGetBound (bound,Eerg)
     EeV = Eerg * op_erg2eV

     if (EeV < zero) then
         call Driver_abortFlash ('[op_setEnergyGroupBoundaries] ERROR: Energy boundary < 0')
     end if

     op_energyGroupBoundaries (bound) = EeV
  end do
!
!
!    ...Ready!
!
!
  return
end subroutine op_setEnergyGroupBoundaries
