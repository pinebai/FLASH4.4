      subroutine amr_guardcell(mype,iopt,nlayers,idir,maxNodetype_gcWanted)


! $RCSfile: amr_guardcell.F90,v $
! $Revision: 1.2 $
! $Date: 2004/10/23 22:39:50 $



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information.
! It does this in a couple of steps:
! (1) by restricting data from leaf blocks to their parents
! (2) then by exchanging data between blocks at the same refinement 
!     level
! (3) and then finally by communicating guard cell data from all 
!     coarser blocks to any of their finer neighbors.
!
! Note, the guardcells lying along block edges are guaranteed
! to be correct. This is necessary because guardcell values in
! children bordering coarser neighbors use their parents edge
! guard cells during the interpolation near the edge of that face.
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------
!
! Arguments:
!      mype local processor number
!
!------------------------------------


      use physicaldata
      use Driver_interface, ONLY : Driver_abortFlash
      use paramesh_interfaces, ONLY : amr_restrict

      use workspace
      use tree
      implicit none
      include 'mpif.h'


#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in) :: mype,iopt,nlayers,idir
      integer,OPTIONAL,intent(IN) :: maxNodetype_gcWanted
      integer :: maxNodetype_gcWanted_loc

      logical ltype2only
      logical :: first = .TRUE.

      integer idiag,ierr, errorcode

      real tot_time
      integer iempty,nprocs




      save first

!------------------------------------



#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting GUARDCELL '
         print *,' starting GUARDCELL '
         close (30)
      end if
#endif
      idiag = 1

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)

      tot_time = 0.

      if (first) then
         call get_tree_nodetypes(nprocs,mype)
         first = .FALSE.
      end if


      
! make sure that nlayers and iopt are set consistently.
      if(iopt.eq.1.and.nlayers.ne.nguard) then
        if(mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
          write(30,*) 'PARAMESH ERROR !'
          write(30,*) 'Error in guardcell - iopt and nlayers'
          write(30,*) 'are not consistent. For iopt=1 you must'
          write(30,*) 'set nlayers=nguard.'
          close(30)
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt=1 you must'
          write(*,*) 'set nlayers=nguard.'
        endif
        call Driver_abortFlash("Error in guardcell: for iopt=1 you must set nlayers=nguard")
      elseif(iopt.eq.2.and.nlayers.gt.nguard_work) then
        if(mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
          write(30,*) 'PARAMESH ERROR !'
          write(30,*) 'Error in guardcell - iopt and nlayers'
          write(30,*) 'are not consistent. For iopt=2 you must'
          write(30,*) 'set nlayers le nguard_work.'
          close(30)
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt=2 you must'
          write(*,*) 'set nlayers le nguard_work.'
        endif
        call Driver_abortFlash("Error in guardcell: for iopt=2 you must set nlayers le nguard")
     endif

     if (present(maxNodetype_gcWanted)) then
        maxNodetype_gcWanted_loc = maxNodetype_gcWanted
     else
        maxNodetype_gcWanted_loc = 1 !by default, guard cells are wanted for leaf blocks only!?
     end if

      if (idiag.eq.1) then
! put recognizably bad data in corners to permit recognition of
! corner guard cells diagonally opposite coarser blocks.
        call amr_mark_edges(mype,iopt)
      end if

! restrict data from leaf blocks to their parents. This will enable the
! parent blocks to provide guard cell data to their neighbors in cases
! where the leaf blocks have coarser neighbors.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting RESTRICT '
         print *,' starting RESTRICT '
         close (30)
      end if
#endif

      iempty = 0
      call amr_restrict(mype,iopt,iempty)

! Apply boundary conditions by filling guard cells at external boundaries;
! copied up here from after amr_guardcell_srl call.

      if (idiag==1) then
#ifdef DEBUG_AMR
         if (mype.eq.0) then
            open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
            write (30,*) ' starting TOT_BND (before GUARDCELL_SRL)'
            print *,' starting TOT_BND  (before GUARDCELL_SRL)'
            close (30)
         end if
#endif

         if (iopt.eq.1) then
            call tot_bnd(idiag,idir)
         else
            call tot_bnd_work(idiag,idir)
         end if
      end if

! blocks provide guard cell data to all their neighbors which share their
! level of refinement.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting GUARDCELL_SRL '
         print *,' starting GUARDCELL_SRL '
         close (30)
      end if
#endif
      
      call amr_guardcell_srl(mype,iopt,nlayers,idiag,idir,maxNodetype_gcWanted_loc)

! apply boundary conditions by filling guard cells at external boundaries.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting TOT_BND '
         print *,' starting TOT_BND '
         close (30)
      end if
#endif

      ltype2only = .false.
      if (iopt.eq.1) then
        call tot_bnd(idiag,idir)
      else
        call tot_bnd_work(idiag,idir)
      end if

! guard cell data is sent to any neighbor blocks with finer resolution.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' starting GUARDCELL_C_TO_F '
         print *,' starting GUARDCELL_C_TO_F '
         close (30)
      end if
#endif
      call amr_guardcell_c_to_f(mype,iopt,nlayers,idiag,idir)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown', & 
     &        position='append')
         write (30,*) ' done GUARDCELL '
         write (30,*) ' '
         print *,' done GUARDCELL '
         print *,' '
         close (30)
      end if
#endif



      return
    end subroutine amr_guardcell




