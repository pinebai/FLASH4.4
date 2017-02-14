!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_finalize
!!
!! NAME
!!
!!  RadTrans_finalize
!!
!! SYNOPSIS
!!
!!  call RadTrans_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the RadTrans unit.
!!
!! ARGUMENTS
!!
!!***
subroutine RadTrans_finalize ()
  use rt_data, ONLY: rt_mgdBounds, rt_mgdDomainBC, rt_mgdBcVals
  implicit none
  
  if(allocated(rt_mgdBounds))   deallocate(rt_mgdBounds)
  if(allocated(rt_mgdDomainBC)) deallocate(rt_mgdDomainBC)
  if(allocated(rt_mgdBcVals))   deallocate(rt_mgdBcVals)

  return
end subroutine RadTrans_finalize
