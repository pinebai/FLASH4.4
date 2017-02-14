!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_copy_soln(level)



!------------------------------------------------------------------------
!
! This routine copies a global solution update from the time 
! synchronized global solution arrays, into the arrays used
! during the solution update, as is required when 
! using NO_PERMANENT_GUARDCELLS and the amr_1blk_guardcell routines.
!
! Arguments:
!      level        integer           if -1 then blocks at all refinement
!                                     levels are copied, otherwise only blocks
!                                     at level are copied.
!
! Written :     Peter MacNeice          May 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use timings
      use tree
      use paramesh_comm_data, ONLY : amr_mpi_meshComm
      use paramesh_mpi_interfaces, only : & 
     &                       mpi_amr_global_domain_limits, & 
     &                       mpi_amr_morton_limits, & 
     &                       mpi_morton_bnd, & 
     &                       mpi_amr_gsurr_blks

      implicit none

      include 'mpif.h'

      integer, intent(in) :: level

      integer :: mype,nprocs,lb,ierr
      integer :: tag_offset, ivar
      integer :: nguard0
      double precision :: time1

!-------------------------

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

        nguard0 = nguard*npgs

        if (no_permanent_guardcells) then

        Call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)


        if(level.eq.-1) then

        do lb = 1, lnblocks

        if(nvar.gt.0) then
          do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then
          gt_unk(ivar,:,:,:,lb) = unk(ivar,:,:,:,lb)
          endif
          enddo
        endif
#ifndef LIBRARY
#ifdef NO_PERMANENT_GUARDCELLS
        if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
          gt_facevarx(ivar,:,:,:,lb) =  & 
     &       facevarx(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
          gt_facevary(ivar,:,:,:,lb) =  & 
     &       facevary(ivar,:,:,:,lb)
          endif
          enddo
          end if
          if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
          gt_facevarz(ivar,:,:,:,lb) =  & 
     &       facevarz(ivar,:,:,:,lb)
          endif
          enddo
          end if
        endif
#endif
#else /*LIBRARY*/
        if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
          gt_facevarx(ivar,:,:,:,lb) =  & 
     &       facevarx(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
          gt_facevary(ivar,:,:,:,lb) =  & 
     &       facevary(ivar,:,:,:,lb)
          endif
          enddo
          end if
          if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
          gt_facevarz(ivar,:,:,:,lb) =  & 
     &       facevarz(ivar,:,:,:,lb)
          endif
          enddo
          end if
        endif
#endif /*LIBRARY*/
        if(nvaredge.gt.0) then
          if (ndim > 1) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(1,ivar)) then
          gt_unk_e_x(ivar,:,:,:,lb) =  & 
     &       unk_e_x(ivar,:,:,:,lb)
          endif
          enddo
          do ivar=1,nvaredge
          if(int_gcell_on_ec(2,ivar)) then
          gt_unk_e_y(ivar,:,:,:,lb) =  & 
     &       unk_e_y(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim == 3) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(3,ivar)) then
          gt_unk_e_z(ivar,:,:,:,lb) =  & 
     &       unk_e_z(ivar,:,:,:,lb)
          endif
          enddo
          end if
          end if
        endif

        if(nvarcorn.gt.0) then
          do ivar=1,nvarcorn
          if(int_gcell_on_nc(ivar)) then
          gt_unk_n(ivar,:,:,:,lb) = unk_n(ivar,:,:,:,lb)
          endif
          enddo
        endif

        end do ! end loop over blocks

      else

        do lb=1,lnblocks
        if(lrefine(lb).eq.level) then

          if(nvar.gt.0) then
          do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then
            gt_unk(ivar,:,:,:,lb) = unk(ivar,:,:,:,lb)
          endif
          enddo
          endif
#ifndef LIBRARY
#ifdef NO_PERMANENT_GUARDCELLS
          if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
            gt_facevarx(ivar,:,:,:,lb) = facevarx(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
            gt_facevary(ivar,:,:,:,lb) = facevary(ivar,:,:,:,lb)
          endif
          enddo
          end if
          if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
            gt_facevarz(ivar,:,:,:,lb) = facevarz(ivar,:,:,:,lb)
          endif
          enddo
          end if
          endif
#endif
#else /*LIBRARY*/
          if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
            gt_facevarx(ivar,:,:,:,lb) = facevarx(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
            gt_facevary(ivar,:,:,:,lb) = facevary(ivar,:,:,:,lb)
          endif
          enddo
          end if
          if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
            gt_facevarz(ivar,:,:,:,lb) = facevarz(ivar,:,:,:,lb)
          endif
          enddo
          end if
          endif
#endif
          if(nvaredge.gt.0) then
          if (ndim > 1) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(1,ivar)) then
            gt_unk_e_x(ivar,:,:,:,lb) = unk_e_x(ivar,:,:,:,lb)
          endif
          enddo
          do ivar=1,nvaredge
          if(int_gcell_on_ec(2,ivar)) then
            gt_unk_e_y(ivar,:,:,:,lb) = unk_e_y(ivar,:,:,:,lb)
          endif
          enddo
          if (ndim == 3) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(3,ivar)) then
            gt_unk_e_z(ivar,:,:,:,lb) = unk_e_z(ivar,:,:,:,lb)
          endif
          enddo
          end if
          end if
          endif
          if(nvarcorn.gt.0) then
          do ivar=1,nvarcorn
          if(int_gcell_on_nc(ivar)) then
            gt_unk_n(ivar,:,:,:,lb) = unk_n(ivar,:,:,:,lb)
          endif
          enddo
          endif

        endif
        enddo

      endif                           ! end of level iftest

      Call MPI_BARRIER(amr_mpi_meshComm, ierr)

! Since this routine is a prelude to calling amr_1blk_guardcell
! we check here to make sure that surrblks has been computed and
! stored for each block.
      if(gsurrblks_set.ne.1) then
        Call MPI_COMM_RANK(amr_mpi_meshComm, mype, ierr)

!--------
! This is a temporary section.
! call mpi_amr_gsurr_blks immediately after call to mpi_morton_bnd
! because it uses pe_source and r_mortonbnd which are reset in the
! other morton_bnd_?? routines. Note this is temporary. Will
! redesign so mpi_amr_gsurr_blks is less context sensitive.

! Find the coordinate ranges
         call mpi_amr_global_domain_limits
!
! Compute and save morton number range for each processor
      call mpi_amr_morton_limits(mype)
!
! Set up surrounding blocks of all local blocks (must not precede
! setting of grid_xmin,... etc)
      tag_offset = 100
      if(nprocs.gt.1) & 
     &    call mpi_morton_bnd(mype,nprocs,tag_offset)
!--------
        Call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)
        call mpi_amr_gsurr_blks(mype,nprocs)

      endif

      else  ! no_permanent_guardcells

      if (force_consistency) then
      if(nfacevar.gt.0) then
        do lb = 1,lnblocks
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
          gt_facevarx(ivar,1,:,:,lb) = facevarx(ivar,1+nguard0,:,:,lb)
          gt_facevarx(ivar,2,:,:,lb) =  & 
     &                             facevarx(ivar,nxb+1+nguard0,:,:,lb)
          endif
          enddo

          if(ndim.ge.2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
          gt_facevary(ivar,:,1,:,lb) =  & 
     &                             facevary(ivar,:,1+nguard0*k2d,:,lb)
          gt_facevary(ivar,:,1+k2d,:,lb) = & 
     &                       facevary(ivar,:,nyb+(1+nguard0)*k2d,:,lb)
          endif
          enddo
          endif

          if(ndim.eq.3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
          gt_facevarz(ivar,:,:,1,lb) =  & 
     &                           facevarz(ivar,:,:,1+nguard0*k3d,lb)
          gt_facevarz(ivar,:,:,1+k3d,lb) = & 
     &                     facevarz(ivar,:,:,nzb+(1+nguard0)*k3d,lb)
          endif
          enddo
          endif

        enddo
      endif
      endif

      endif ! no_permanent_guardcells

      Call MPI_BARRIER(amr_mpi_meshComm, ierr)

      if (timing_mpi) then
              timer_amr_1blk_copy_soln =  timer_amr_1blk_copy_soln & 
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine amr_1blk_copy_soln
