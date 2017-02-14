!!****if* source/Grid/GridParticles/GridParticlesMove/Sieve/BlockMatch/gr_ptResetProcPair
!!
!! NAME
!!  gr_ptResetProcPair
!!
!! SYNOPSIS
!!
!!  gr_ptResetProcPair( logical(IN)  :: haveNonLocalData
!!                      logical(OUT) :: mustCommunicate)
!!
!!  
!! DESCRIPTION 
!!  
!!  This routine resets the sieve iterations and finds out if any parallel
!!   processing of the sieve is needed at all
!!
!! ARGUMENTS 
!!
!!   haveNonLocalData : indication that there is some non local data present at myPE
!!   mustCommunicate : indicates that the iterations of seive movement must continue
!!                     because some processor has non local data
!!
!!
!! NOTES
!!   
!!
!!
!!
!!***

subroutine gr_ptResetProcPair(haveNonLocalData, mustCommunicate)
  use Driver_interface, ONLY : Driver_abortFlash

  use gr_ptData, ONLY : gr_ptToggle, gr_ptProcDist, gr_ptSieveFreq, gr_ptSieveCheckFreq


  use Grid_data, ONLY : gr_meshNumProcs

  implicit none

  include "Flash_mpi.h"

  logical, intent(IN)  :: haveNonLocalData
  logical, intent(OUT) :: mustCommunicate
  
  integer :: ierr

  gr_ptProcDist = 1
  gr_ptToggle = 1
  gr_ptSieveCheckFreq=gr_ptSieveFreq+1
  
  mustCommunicate = .false.
  if(gr_meshNumProcs==1) then
     write(*,*) 'haveNonLocalData',haveNonLocalData,' , mustCommunicate=',mustCommunicate
     if(haveNonLocalData) call Driver_abortFlash("gr_sieve: destCount>0 for 1 proc")
  else
     call MPI_ALLREDUCE(haveNonLocalData, mustCommunicate, 1, FLASH_LOGICAL, &
             FLASH_LOR, FLASH_COMM, ierr)  
  end if

end subroutine gr_ptResetProcPair

