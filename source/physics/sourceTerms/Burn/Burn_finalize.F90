!!****f* source/physics/sourceTerms/Burn/Burn_finalize
!!
!! NAME
!!  
!!  Burn_finalize
!!
!!
!! SYNOPSIS
!! 
!!  call Burn_finalize()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the Burn module.
!!
!! NOTES
!!  
!!  There is no implementation that does anything.
!!
!!  The NSE arrays used with parametricBurn are deallocated by
!!  NSE_finalize, which should be called directly frm Driver_finalizeFlash.
!!
!!***


subroutine Burn_finalize()


  implicit none

  return

end subroutine Burn_finalize
