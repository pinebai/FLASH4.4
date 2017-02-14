!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_checkGridState
!!
!! NAME
!!  gr_checkGridState
!!
!! SYNOPSIS
!!
!!  gr_checkGridState()
!!
!! DESCRIPTION
!!
!!  Checks that the cached grid data contains valid values.
!!  If we find invalid (block,proc) identifiers then we declare the 
!!  grid to be in an inconsistent state, and we abort.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!! 
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_checkGridState()

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interfaceTypeDecl
  use gr_interface, ONLY : gr_findAllNeghID
  implicit none

  include "Flash_mpi.h"
  integer, dimension(MAXBLOCKS) :: listofBlocks
  integer :: blkCount, blk, blockID, allCenters
  integer :: ib, ie, jb, je, kb, ke, i, j, k, numNegh, ierr, eachNegh
  integer, dimension(LOW:HIGH,MDIM) :: guardCoords
  integer, dimension(MDIM) :: guardCellID, blkSize, guard
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  logical :: allNeghValid
  type (AllBlockRegions_t) :: surrBlksSummary
  integer, dimension(BLKNO:PROCNO) :: neghBlkProc

  allNeghValid = .true.
  blockID = 1
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1
  guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)


  call Grid_getListOfBlocks(LEAF, listofBlocks, blkCount)
  do blk = 1, blkCount
     blockID = listofBlocks(blk)

     call gr_findAllNeghID(blockID, surrBlksSummary)

     !Loop over each guard cell region.          
     allCenters = 2**NDIM
     kb=blkLimitsGC(LOW,KAXIS); ke=kb+K3D*(guard(KAXIS)-1)
     do k=LEFT_EDGE,LEFT_EDGE+K3D*(RIGHT_EDGE-1)
        jb=blkLimitsGC(LOW,JAXIS); je=jb+K2D*(guard(JAXIS)-1)
        do j=LEFT_EDGE,LEFT_EDGE+K2D*(RIGHT_EDGE-1)
           ib=blkLimitsGC(LOW,IAXIS); ie=ib+guard(IAXIS)-1
           do i=LEFT_EDGE,RIGHT_EDGE
              if((i*j*k)/=allCenters) then !! make sure it is gCell region

                 guardCellID(IAXIS)=i; guardCellID(JAXIS)=j; guardCellID(KAXIS)=k
                 guardCoords(LOW,IAXIS)=ib; guardCoords(HIGH,IAXIS)=ie
                 guardCoords(LOW,JAXIS)=jb; guardCoords(HIGH,JAXIS)=je
                 guardCoords(LOW,KAXIS)=kb; guardCoords(HIGH,KAXIS)=ke

                 !IN: blockID, guardCellID
                 !OUT: arrNeghID,neghCornerID,numNegh
                 !call gr_ptFindNegh(blockID,guardCellID,arrNeghID,neghCornerID,numNegh)
                 numNegh = surrBlksSummary % regionInfo(i,j,k) % numNegh


                 !Ensures we handle numNegh=0, i.e. non-periodic global boundary correctly.
                 if (numNegh > 0) then
                    do eachNegh = 1, numNegh

                       neghBlkProc(BLKNO:PROCNO) = surrBlksSummary &
                            % regionInfo(i,j,k) % details(BLKNO:PROCNO,eachNegh)

                       if ( (neghBlkProc(BLKNO) <= 0) .or. &
                            (neghBlkProc(BLKNO) > MAXBLOCKS) .or. &
                            (neghBlkProc(PROCNO) < 0) .or. & 
                            (neghBlkProc(PROCNO) >= gr_meshNumProcs) ) then
                          
                          allNeghValid = .false.
                          print *, "Neighbor to block:", blockID, &
                               "on processor:", gr_meshMe, &
                               "is invalid.  Details:", neghBlkProc
                          
                       end if
                    end do                  
                 end if

              end if !End of if all centers
              ib=ie+1;ie=ib+mod(i,CENTER)*blkSize(IAXIS)+mod(i+1,CENTER)*guard(IAXIS)-1
           end do !End of do i
           jb=je+1;je=jb+mod(j,CENTER)*blkSize(JAXIS)+mod(j+1,CENTER)*guard(JAXIS)-1
        end do !End of do j
        kb=ke+1;ke=kb+mod(k,CENTER)*blkSize(KAXIS)+mod(k+1,CENTER)*guard(KAXIS)-1
     end do !End of do k

  end do  


  call MPI_Barrier(gr_meshComm, ierr)
  if (allNeghValid .eqv. .false.) then
     call Driver_abortFlash("Invalid neighbors found!")
  end if

  !So we only reach success message when all processors are valid.
  call MPI_Barrier(gr_meshComm, ierr)  
  if (gr_meshMe == 0) then
     print *, "All neighbors valid"
  end if

end subroutine gr_checkGridState  
