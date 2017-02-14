!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_finalize
!!
!! NAME
!!
!!  TreeRay_finalize
!!  
!! SYNOPSIS
!!
!!  TreeRay_finalize()
!!
!! DESCRIPTION
!!  deallocate any memory that might have allocated in TreeRay unit
!!  and other housekeeping to prepare for the end of the unit.
!!
!!
!!***
subroutine TreeRay_finalize()
  use TreeRay_data, ONLY : tr_bhMassRays, tr_bhVolRays, tr_bhSrcfRays, &
    tr_intersectList, tr_bhCdMaps, tr_bhUseTreeRay
  implicit none

  if (.not. tr_bhUseTreeRay) return

  deallocate(tr_intersectList)

#ifdef TR_OPTICALDEPTH
  deallocate(tr_bhCdMaps)
#endif


  return

end subroutine TreeRay_finalize
