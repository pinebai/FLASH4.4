!!****if* source/Grid/GridParticles/GridParticlesMove/Sieve/BlockMatch/gr_ptNextProcPair
!!
!! NAME
!!  gr_ptNextProcPair
!!
!! SYNOPSIS
!!
!!  gr_ptNextProcPair(  integer(IN)  :: timesInLoop,
!!                      integer(OUT) :: src, 
!!                      integer(OUT) :: dest,
!!                      logical(IN)  :: haveNonLocalData
!!                      logical(OUT) :: mustCommunicate)
!!
!!  
!! DESCRIPTION 
!!  
!!  This routine finds the pair of processors with which MyPe must communicate
!!  in the next iteration of sieve motion. The routine is abstrated out from 
!!  gr_ptMoveSieve to allow different algorithms for picking the next processor pair.
!!  For example in some situations the processor pair may remain fixed, as when sorting
!!  processors based on their tags among processors. In the default mode, the number the 
!!  motion is back-and-forth with increasing intervals
!!
!! ARGUMENTS 
!!
!!   timesInLoop : the count of iterations of sieve movement
!!   src : the negh to receive the sieve from
!!   dest : the nehg to send the sieve to
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

#define DEBUG_GRIDNEDATA

subroutine gr_ptNextProcPair(timesInLoop, src, dest, haveNonLocalData, mustCommunicate)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use gr_ptData, ONLY : gr_ptToggle, gr_ptProcDist, &
       gr_ptSieveFreq, gr_ptSieveCheckFreq
       

  implicit none

  include "Flash_mpi.h"  
  integer, intent(IN)  :: timesInLoop
  integer, intent(OUT) :: src, dest
  logical, intent(IN)  :: haveNonLocalData
  logical, intent(OUT) :: mustCommunicate
  integer :: ierr

  if(timesInLoop > gr_meshNumProcs) then
     if(haveNonLocalData) then
        print*,'I still have data for communication'
     end if
     call Driver_abortFlash("gr_sieve: timesInLoop > gr_meshNumProcs, neData should have been filtered by now")
  else
     gr_ptSieveCheckFreq=gr_ptSieveCheckFreq-1
     if((gr_ptSieveCheckFreq==0).or.(timesInLoop==gr_meshNumProcs)) then
        gr_ptSieveCheckFreq=gr_ptSieveFreq
        mustCommunicate = .false.
        call MPI_ALLREDUCE(haveNonLocalData, mustCommunicate, 1, MPI_LOGICAL, &
             MPI_LOR, gr_meshComm, ierr)  
     else
        mustCommunicate=.true.
     end if
  end if

  if(mustCommunicate) then
     dest = gr_meshMe + gr_ptProcDist*gr_ptToggle
     src = gr_meshMe - gr_ptProcDist*gr_ptToggle
     
     
     !catch ending conditions
     if(dest >= gr_meshNumProcs) then
        dest = dest - gr_meshNumProcs
     end if
     
     if(src < 0) then
        src = src + gr_meshNumProcs
     end if
     
     if(src >= gr_meshNumProcs) then
        src = src - gr_meshNumProcs
     end if
     
     if(dest < 0) then
        dest = dest + gr_meshNumProcs
     end if

     
#ifdef DEBUG_GRIDNEDATA
     if(dest < 0) then
        call Driver_abortFlash("Grid_moveNeData: dest < 0")
     end if
     
     if(dest >= gr_meshNumProcs) then
        call Driver_abortFlash("Grid_moveNeData: dest >= gr_meshNumProcs")
     end if
     
     if(src < 0) then
        call Driver_abortFlash("Grid_moveNeData: src < 0")
     end if
     
     if(src >= gr_meshNumProcs) then
        call Driver_abortFlash("Grid_moveNeData: scr >= gr_meshNumProcs")
     end if
     
     if(src == gr_meshMe) then
        call Driver_abortFlash("src == gr_meshMe, and we still have unmatched neData, bad")
     end if
     
     if(dest == gr_meshMe) then
        call Driver_abortFlash("dest == gr_meshMe, and we still have unmatched neData, bad")
     end if
     
#endif
     
     !set dest and src for next time through loop

     gr_ptToggle = gr_ptToggle * (-1)     
     gr_ptProcDist = gr_ptProcDist + 1

  else
     src=gr_meshMe
     dest=gr_meshMe
  end if
  
end subroutine gr_ptNextProcPair
