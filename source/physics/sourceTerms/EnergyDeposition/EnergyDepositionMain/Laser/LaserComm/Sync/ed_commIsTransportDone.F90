!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Sync/ed_commIsTransportDone
!!
!!  NAME     
!!   ed_commIsTransportDone
!!
!!  SYNOPSIS
!!   ed_commIsTransportDone(logical, intent(OUT) :: isTransportDone)
!!
!!  DESCRIPTION 
!!
!!  ARGUMENTS
!!
!!  NOTES
!!
!!  SIDE EFFECTS
!!
!!***

#include "constants.h"
#include "EnergyDeposition.h"

subroutine ed_commIsTransportDone(isTransportDone)
  use EnergyDeposition_data, ONLY : ed_rayCount, ed_meshComm, ed_rays
  implicit none
  include "Flash_mpi.h"
  logical, intent(OUT) :: isTransportDone
  integer :: error
  logical :: activeRays, moreRays

  moreRays = (ed_rayCount > 0)
  if (moreRays) then
     moreRays = any (ed_rays (RAY_BLCK,1:ed_rayCount) /= real (NONEXISTENT))
  end if
  call MPI_Allreduce (moreRays,      &
       activeRays,    &
       1,             &
       MPI_LOGICAL,   &
       MPI_LOR,       &
       ed_meshComm,   &
       error          )
  isTransportDone = .not.activeRays
end subroutine ed_commIsTransportDone
