!!****f* source/Driver/Driver_superTimeStep
!!
!! NAME
!!
!!  Driver_superTimeStep
!!
!! SYNOPSIS
!!
!!  Driver_superTimeStep()
!!
!! DESCRIPTION
!!
!! This routine implements the super time steppping advancement algorithm
!! to overcome small diffusive time scales in explicit formulation
!!
!!  
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_superTimeStep(dt,nuSTS,nstepSTS,nstepTotalSTS,dt_subSTS)

  implicit none
  real, intent(IN)    :: dt,nuSTS
  integer, intent(IN) :: nstepSTS,nstepTotalSTS
  real, intent(OUT)   :: dt_subSTS

  return
end subroutine Driver_superTimeStep
