!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_checkGridConsistency
!!
!! NAME
!!  gr_checkGridConsistency
!!
!! SYNOPSIS
!!
!!  call gr_checkGridConsistency()
!!
!! DESCRIPTION
!!
!!  Checks that the grid is consistent.
!!  If we find invalid refinement level changes, then we declare the 
!!  grid to be in an inconsistent state, and we abort.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!! 
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_checkGridConsistency()
  use IO_interface, ONLY : IO_writePlotfile
  use Driver_interface, ONLY: Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stamp
  use paramesh_mpi_interfaces, ONLY: mpi_amr_local_surr_blks_lkup
  use tree, ONLY : nodetype, lnblocks
  use Grid_data, ONLY: gr_meshMe, gr_meshNumProcs, gr_meshComm
  implicit none
#include "Flash_mpi.h"
  logical, parameter :: mayAttemptFix = .TRUE.
  logical :: attemptedFix
  logical :: missingNeighbor, linconsistent, ginconsistent
  integer :: lb, ierr, i,j,k
  Integer :: surrblks(3,3,3,3)
  Integer :: psurrblks(3,3,3,3)

  attemptedFix = .FALSE.
1 CONTINUE
  linconsistent = .FALSE.
  do lb = 1, lnblocks
     call mpi_amr_local_surr_blks_lkup(gr_meshMe,lb,                       & 
                                surrblks,.FALSE.,psurrblks)
     missingNeighbor = .FALSE.
     do k = 2-K3D,2+K3D
        do j = 2-K2D,2+K2D
           do i = 1,3

              if (surrblks(1,i,j,k) > -20.AND.surrblks(1,i,j,k) < 0)       & 
                   missingNeighbor = .TRUE.
           end do
        end do
     end do
     if (missingNeighbor .AND. nodetype(lb) .NE. 1) then
        linconsistent = .TRUE.
#if NDIM == 1
#define THREE_POW 3
#elif NDIM == 2
#define THREE_POW 9
#else
#define THREE_POW 27
#endif
999     format("Missing neighbors on PE",I5,", lb=",I5,",surrblks=", &
             THREE_POW("(",I4,",",I4,",",I2,")"))
        print 999,gr_meshMe,lb,surrblks(1:3,:,2-K2D:2+K2D,2-K3D:2+K3D)
     end if
  end do

  call MPI_ALLREDUCE(linconsistent, ginconsistent, 1, MPI_LOGICAL, MPI_LOR, gr_meshComm, ierr)
  if  (ginconsistent .AND. mayAttemptFix .AND. .NOT. attemptedFix) then
     call Logfile_stamp(&
          "Detected inconsistent Grid state, FLASH will attempt fix...",&
          "[gr_checkGridConsistency]")
     call gr_ensureValidNeighborInfo(0)
     attemptedFix = .TRUE.
     goto 1
  end if
     
  if  (ginconsistent) then
     call Logfile_stamp(&
          "Detected inconsistent Grid state, FLASH will dump a plot file for analysis and then abort!",&
          "[gr_checkGridConsistency]")
     call IO_writePlotfile( .true.)     
     call MPI_BARRIER(gr_meshComm, ierr)
     call dr_sleep(5)
     if (linconsistent) call Driver_abortFlash("missing neighbor!")
  end if

end subroutine gr_checkGridConsistency
