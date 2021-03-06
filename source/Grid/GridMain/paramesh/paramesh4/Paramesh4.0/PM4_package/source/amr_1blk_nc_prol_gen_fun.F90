!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_nc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype)


!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from one of the arrays unk_n1, and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written directly into unk_n1.
! The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! This particular prolongation is simple linear interpolation. It can
! only be used for blocks with an even number of grid cells.
!
! Conservative prolongation. Special treatment for the  cells immediately
! adjacent to a boundary (ie i=nguard,nguard+1,iu_bnd1-nguard,iu_bnd1-nguard+1
! and likewise for j and k indeces) if using an even number of grid cells
! per block along that axis. No special treatment is required when the number
! of cells is odd.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_face_fun_init has been called.
!
!
! Written :     Peter MacNeice          December 2000
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use prolong_arrays

      use paramesh_interfaces, only : amr_1blk_nc_prol_linear, & 
     &                                amr_1blk_nc_prol_genorder, & 
     &                                amr_1blk_nc_prol_user


      implicit none

!------------------------------------

      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:,:)

      integer :: ivar


      include 'mpif.h'
      double precision :: time1
      
!------------------------------------

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      do ivar=1,nvarcorn

      if (int_gcell_on_nc(ivar)) then

      if (interp_mask_nc(ivar) < 20) then

      if (interp_mask_nc(ivar) == 1) then
! linear interpolation
         call amr_1blk_nc_prol_linear(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype,ivar)
      else
! general interpolation routine
         call amr_1blk_nc_prol_genorder(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype,ivar,interp_mask_nc(ivar))
      end if

      else
         
         call amr_1blk_nc_prol_user()

      end if

      end if

      end do

      if (timing_mpi) then
              timer_amr_1blk_nc_prol_gen =  & 
     &                          timer_amr_1blk_nc_prol_gen & 
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine amr_1blk_nc_prol_gen_fun





