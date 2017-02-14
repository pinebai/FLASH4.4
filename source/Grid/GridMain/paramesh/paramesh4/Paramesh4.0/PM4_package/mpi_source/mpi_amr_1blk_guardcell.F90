!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_1blk_guardcell
!! NAME
!!
!!   amr_1blk_guardcell
!!
!! SYNOPSIS
!!
!!   call amr_1blk_guarcell(mype, iopt, nlayers, lb, pe, 
!!                          lcc, lfc, lec, lnc, l_srl_only,
!!                          icoord, ldiag)
!!   call amr_1blk_guarcell(mype, iopt, nlayers, lb, pe, 
!!                          lcc, lfc, lec, lnc, l_srl_only,
!!                          icoord, ldiag, 
!!                          nlayersx, nlayersy, nlayersz)
!!
!!   call amr_1blk_guarcell(integer, integer, integer, integer, integer, 
!!                          logical, logical, logical, logical, logical,
!!                          integer, logical, 
!!                          optional integer, optional integer, optional integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: iopt           
!!        A switch to control which data source is to be used:
!!         iopt=1 will use 'unk', 'facevarx', 'facevary', 'facevarz', 
!!                'unk_e_x', 'unk_e_y', 'unk_e_z', and 'unk_n'
!!         iopt>=2 will use 'work'
!!
!!   integer, intent(in) :: nlayers        
!!        The number of guard cell layers at each boundary 
!!         (Note: this argument is deprecated not longer functions !).
!!
!!   integer, intent(in) :: lb             
!!          The block number on processor pe selected for guardcell filling.
!!
!!   integer, intent(in) :: pe             
!!        The processor storing the block selected for guarcell filling.
!!
!!   logical, intent(in) :: lcc            
!!        A logical switch controlling whether unk or work data is filled.
!!
!!   logical, intent(in) :: lfc            
!!        A logical switch controlling whether facevarx(y)(z) data is filled.
!!
!!   logical, intent(in) :: lec            
!!        A logical switch controlling whether unk_e_x(y)(z) data is filled.
!!
!!   logical, intent(in) :: lnc            
!!        A logical switch controlling whether unk_n data is filled.
!!
!!   logical, intent(in) :: l_srl_only     
!!        A logical switch which, if true, switches off the filling from 
!!        coarse neighbors. This is used during restriction when odd  
!!        block sizes are in use.
!!
!!   integer, intent(in) :: icoord         
!!        An integer switch used to select which faces of the block are to be 
!!        considered. If icoord=0 all faces are considered. If icoord=1 only 
!!        faces perp. to the x-axis are considered, if icoord=2 only faces
!!        perp. to the y-axis are considered, and if icoord=3 only faces perp. 
!!        to the z-axis are considered.
!!
!!   logical, intent(in) :: ldiag          
!!        A logical switch which controls whether guardcells corresponding to 
!!        neighbor blocks diagonally opposite block edges and corners are filled.
!!
!!   optional, integer, intent(in) :: nlayersx, nlayersy, nlayersz
!!        Optional arguments to specify the number of guardcells to fill in each 
!!        coordninate direction.  If not specified, then the default is to fill all
!!        available guardcells.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   workspace
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!   paramesh_interfaces
!!
!! CALLS
!!
!!     mpi_amr_local_surr_blks_lkup
!!     mpi_amr_1blk_guardcell_c_to_f
!!     mpi_amr_get_remote_block
!!     amr_perm_to_1blk
!!     amr_1blk_guardcell_srl
!!     amr_1blk_guardcell_f_to_c
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells are filled with data for 
!!   the block secified in the call to this routine.
!!
!! DESCRIPTION
!!
!!   This routine manages the transfer of guard cell data for a
!!   specific single block. It uses the morton numbering scheme for the
!!   determination of the neighbor relations and differs in that from
!!   the older routine amr_1blk_guardcell.  The user can choose how many layers 
!!   of guardcells are filled on each face of a block by specifying the optional 
!!   arguments 'nlayersx', 'nlayersy', and 'nlayersz'.
!!
!! NOTES
!!
!!  This routine must be used with care !
!!  This routine was written to be used in a code as illustrated
!!  in the following snippet of pseudo-code
!!
!!              .
!!              .
!!              .
!!        synchronization pt
!!
!!        loop over grid blocks
!!          (copy current block into working block and fill its guardcells)
!!          (perform some set of operations on block)
!!          (store result from working block)
!!        end loop
!!
!!        synchronization pt
!!              .
!!              .
!!              .
!!
!!  Caveat 1:
!!   If you are using this routine, you must remember to save the solution
!!   at the first synchronization point (ie call amr_1blk_copy_soln), so 
!!   that each block uses the same time synchronized solution during its 
!!   guardcell filling.
!!
!!  Caveat 2:
!!   It is implicitly assumed that the parent blocks on all leaf nodes
!!   have valid data. (This is necessary to ensure that a general restriction
!!   operator can be supported in the neighborhood of a jump in refinement.)
!!   If ADVANCE_ALL_LEVELS is defined, then this will generally be true. 
!!   However, if the solution is being time-advanced on leaf blocks only 
!!   this may not be true. In this case you should call amr_restrict  
!!   before the first synchronization pt in the example above.
!!   If you are using blocks with an even number of grid points and the
!!   default restriction operators this is not necessary.
!!
!! AUTHORS
!!
!!   Michael Gehmeyr & Peter MacNeice (November 1999) with modifications by 
!!   Kevin Olson for layered guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG
!#define DEBUGX

      subroutine amr_1blk_guardcell( & 
     &                              mype,iopt,nlayers,lb,pe, & 
     &                              lcc,lfc,lec,lnc, & 
     &                              l_srl_only,icoord,ldiag, & 
     &                              nlayersx,nlayersy,nlayersz)




      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace
      use mpi_morton
      use constants
      use paramesh_comm_data

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_amr_local_surr_blks_lkup, & 
     &                                mpi_amr_1blk_guardcell_c_to_f, & 
     &                                mpi_amr_get_remote_block

      use paramesh_interfaces, only : amr_perm_to_1blk, & 
     &                                amr_1blk_guardcell_srl, & 
     &                                amr_1blk_guardcell_f_to_c

      use gr_flashHook_interfaces

      implicit none

      include 'mpif.h'

      logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer, intent(in) :: mype,lb,pe,iopt,nlayers,icoord
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!------------------------------------

! local arrays and variables

      integer ::  nprocs

      integer :: parent_lb,parent_pe
      integer :: i,j,k, ll, idest, iblock
      integer :: ij,ijk

      logical :: lcoarse,l_parent
      logical :: loc_lcc,loc_lfc,loc_lec,loc_lnc
      logical :: lfound

      integer :: surrblks(3,3,3,3), tsurrblks(3,3,3,3)
      integer :: psurrblks(3,3,3,3)
      integer :: pcache_pe,pcache_blk
      integer :: ierrorcode,ierr

      integer :: ipolar(2)
      integer :: ippolar(2)
      integer :: icoord_loc
      logical :: ldiag_loc

      integer :: nlayers0x, nlayers0y, nlayers0z

      logical, save :: first_cc = .true.
      logical, save :: first_nc = .true.
      logical, save :: first_ec = .true.
      logical, save :: first_fc = .true.

      real    :: eps,accuracy,pbnd_box(2,3)

      double precision :: time1
      double precision :: time2

!------------------------------------
!
! The sequence of operations required to fill the guardcells of an
! individual block is :
!
! Step 1:
! Construct a list of blocks surrounding block lb
!
! Step 2:
! Check for coarse neighbors
!
! Step 3:
! Put leaf block data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding.
!
! Step 4:
! Put data from leaf blocks parent into the data_1blk.fh datastructures, 
! with the appropriate guardcell padding. Check to see if this data is 
! currently cached.
!
! Step 5:
! Construct a list of blocks surrounding block lb's parent
!
! Step 6:
! Do guardcell filling for lb's parent from any surrounding blocks at 
! the same refinement level as this parent.
!
! Step 7:
! Do guardcell filling from coarse neigbors into the current block
!
! Step 8:
! Do guardcell filling from any surrounding blocks at the same refinement
! level as the leaf block.
!
! Step 9:
! Apply boundary conditions.
!
!
!------------------------------------

#ifdef DEBUG_FLOW_TRACE
      write(*,*)'amr_1blk_guardcell: pe ',mype,' entered',' iopt ',iopt
#endif /* DEBUG_FLOW_TRACE */

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

      if (timing_mpi) then
         time1 = mpi_wtime()
         time2 = mpi_wtime()
      endif


      call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)

      if (iopt == 1) then
         if (.not.present(nlayersx)) then
            nlayers0x = nguard
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard
         else
            nlayers0z = nlayersz
         end if
      else
#ifdef DEBUG
      write(*,*) 'amr_1blk_guardcell: pe ',mype,' start layers '
#endif /* DEBUG */
         if (.not.present(nlayersx)) then
            nlayers0x = nguard_work
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard_work
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard_work
         else
            nlayers0z = nlayersz
         end if
      end if

#ifdef DEBUG
      write(*,*) 'amr_1blk_guardcell: pe ',mype,' layers1 ended'
#endif /* DEBUG */
      if (nxb/nguard < 2) nlayers0x = min(nlayers0x+1,   nguard)
      if (nyb/nguard < 2) nlayers0y = min(nlayers0y+k2d, nguard)
      if (nzb/nguard < 2) nlayers0z = min(nlayers0z+k3d, nguard)
      
      if (iopt == 1 .and.  & 
     &    all(interp_mask_unk(:) < 20) .and. & 
     &    all(interp_mask_facex(:) < 20) .and. & 
     &    all(interp_mask_facey(:) < 20) .and. & 
     &    all(interp_mask_facez(:) < 20) .and. & 
     &    all(interp_mask_ec(:) < 20) .and. & 
     &    all(interp_mask_nc(:) < 20) ) then
         if (any(nguard < interp_mask_unk(:)-2) .or. & 
     &       any(nguard < interp_mask_facex(:)-2) .or. & 
     &       any(nguard < interp_mask_facey(:)-2) .or. & 
     &       any(nguard < interp_mask_facez(:)-2) .or. & 
     &       any(nguard < interp_mask_ec(:)-2) .or. & 
     &       any(nguard < interp_mask_nc(:)-2)) then
            print *,' PARAMESH ERROR: nguard is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nguard = ',nguard
            print *,' maxval(interp_mask_unk) = ', & 
     &                maxval(interp_mask_unk(:))
            print *,' maxval(interp_mask_facex) = ', & 
     &                maxval(interp_mask_facex(:))
            print *,' maxval(interp_mask_facey) = ', & 
     &                maxval(interp_mask_facey(:))
            print *,' maxval(interp_mask_facez) = ', & 
     &                maxval(interp_mask_facez(:))
            print *,' maxval(interp_mask_ec) = ', & 
     &                maxval(interp_mask_ec(:))
            print *,' maxval(interp_mask_nc) = ', & 
     &                maxval(interp_mask_nc(:))
            call amr_abort()
         end if
         if ( (any(nlayers0x < interp_mask_unk(:)-2) .and. & 
     &             nvar > 0) .or. & 
     &        (any(nlayers0x < interp_mask_facex(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0x < interp_mask_facey(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0x < interp_mask_facez(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0x < interp_mask_ec(:)-2) .and.  & 
     &             nvaredge > 0) .or. & 
     &        (any(nlayers0x < interp_mask_nc(:)-2) .and.  & 
     &             nvarcorn > 0) ) then
            print *,' PARAMESH ERROR: nlayersx is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersx = ',nlayers0x
            print *,' maxval(interp_mask_unk) = ', & 
     &                maxval(interp_mask_unk(:))
            print *,' maxval(interp_mask_facex) = ', & 
     &                maxval(interp_mask_facex(:))
            print *,' maxval(interp_mask_facey) = ', & 
     &                maxval(interp_mask_facey(:))
            print *,' maxval(interp_mask_facez) = ', & 
     &                maxval(interp_mask_facez(:))
            print *,' maxval(interp_mask_ec) = ', & 
     &                maxval(interp_mask_ec(:))
            print *,' maxval(interp_mask_nc) = ', & 
     &                maxval(interp_mask_nc(:))
            call amr_abort()
         end if
         if ( (any(nlayers0y < interp_mask_unk(:)-2) .and.  & 
     &             nvar > 0) .or. & 
     &        (any(nlayers0y < interp_mask_facex(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0y < interp_mask_facey(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0y < interp_mask_facez(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0y < interp_mask_ec(:)-2) .and.  & 
     &             nvaredge > 0) .or. & 
     &        (any(nlayers0y < interp_mask_nc(:)-2) .and.  & 
     &             nvarcorn > 0) ) then
            print *,' PARAMESH ERROR: nlayersy is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersy = ',nlayers0y
            print *,' maxval(interp_mask_unk) = ', & 
     &                maxval(interp_mask_unk(:))
            print *,' maxval(interp_mask_facex) = ', & 
     &                maxval(interp_mask_facex(:))
            print *,' maxval(interp_mask_facey) = ', & 
     &                maxval(interp_mask_facey(:))
            print *,' maxval(interp_mask_facez) = ', & 
     &                maxval(interp_mask_facez(:))
            print *,' maxval(interp_mask_ec) = ', & 
     &                maxval(interp_mask_ec(:))
            print *,' maxval(interp_mask_nc) = ', & 
     &                maxval(interp_mask_nc(:))
            call amr_abort()
         end if
         if ( (any(nlayers0z < interp_mask_unk(:)-2) .and.  & 
     &             nvar > 0) .or. & 
     &        (any(nlayers0z < interp_mask_facex(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0z < interp_mask_facey(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0z < interp_mask_facez(:)-2) .and.  & 
     &             nfacevar > 0) .or. & 
     &        (any(nlayers0z < interp_mask_ec(:)-2) .and.  & 
     &             nvaredge > 0) .or. & 
     &        (any(nlayers0z < interp_mask_nc(:)-2) .and.  & 
     &             nvarcorn > 0) ) then
            print *,' PARAMESH ERROR: nlayersz is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersz = ',nlayers0z
            print *,' maxval(interp_mask_unk) = ', & 
     &                maxval(interp_mask_unk(:))
            print *,' maxval(interp_mask_facex) = ', & 
     &                maxval(interp_mask_facex(:))
            print *,' maxval(interp_mask_facey) = ', & 
     &                maxval(interp_mask_facey(:))
            print *,' maxval(interp_mask_facez) = ', & 
     &                maxval(interp_mask_facez(:))
            print *,' maxval(interp_mask_ec) = ', & 
     &                maxval(interp_mask_ec(:))
            print *,' maxval(interp_mask_nc) = ', & 
     &                maxval(interp_mask_nc(:))
            call amr_abort()
         end if
      elseif (iopt>=2 .and. all(interp_mask_work(:) < 20) ) then
         if (any(nguard_work < interp_mask_work(:)-2)) then
            print *,' PARAMESH ERROR: nguard_work is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nguard_work = ',nguard_work
            print *,' maxval(interp_mask_work) = ', & 
     &                maxval(interp_mask_work(:))
            call amr_abort()
         end if
         if (any(nlayers0x < interp_mask_work(:)-2) .and.  & 
     &           nvar_work > 0) then
            print *,' PARAMESH ERROR: nlayersx is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersx = ',nlayers0x
            print *,' maxval(interp_mask_work) = ', & 
     &                maxval(interp_mask_work(:))
            call amr_abort()
         end if
         if (any(nlayers0y < interp_mask_work(:)-2) .and.  & 
     &           nvar_work > 0) then
            print *,' PARAMESH ERROR: nlayersy is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersy = ',nlayers0y
            print *,' maxval(interp_mask_work) = ', & 
     &                maxval(interp_mask_work(:))
            call amr_abort()
         end if
         if (any(nlayers0z < interp_mask_work(:)-2) .and.  & 
     &           nvar_work > 0) then
            print *,' PARAMESH ERROR: nlayers0z is too small for the  & 
     & chosen interpolation order !!! '
            print *,' nlayersz = ',nlayers0z
            print *,' maxval(interp_mask_work) = ', & 
     &                maxval(interp_mask_work(:))
            call amr_abort()
         end if
      end if
#ifdef DEBUG
      write(*,*) 'amr_1blk_guardcell: pe ',mype,' layers2 ended'
#endif /* DEBUG */

      if(pe.ne.mype) then
          write(*,*) 'Error : trying to fill guardcells for a ', & 
     &               'remote block - not supported with the mpi ', & 
     &               'version. '
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      endif
      if(  & 
     &    mpi_pattern_id.eq.10  .or. & 
     &    (mpi_pattern_id.eq.20.and.lprolong_in_progress) .or. & 
     &    (mpi_pattern_id.eq.40.and.lrestrict_in_progress) & 
     &                            ) then
      elseif(nprocs.gt.1) then
        write(*,*) 'Paramesh error : amr_1blk_guardcell : ', & 
     & ' wrong pattern being', & 
     & ' used for pre-communication for guardcell fill : Fix', & 
     & ' - insert appropriate call to mpi_amr_comm_setup ' & 
     & ,' mpi_pattern_id ',mpi_pattern_id, & 
     & ' lprolong_in_progress ',lprolong_in_progress, & 
     & ' lrestrict_in_progress ',lrestrict_in_progress
        call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      endif

      if (amr_error_checking) then

      if (iopt == 1) then   

        loc_lcc = .false.
        loc_lfc = .false.
        loc_lnc = .false.
        loc_lec = .false.
        if(nvar.gt.0)     loc_lcc = .true.
        if(nfacevar.gt.0) loc_lfc = .true.
        if(nvaredge.gt.0) loc_lec = .true.
        if(nvarcorn.gt.0) loc_lnc = .true.
        if( (lcc.and.(.not.loc_lcc)) .or. (lfc.and.(.not.loc_lfc))  & 
     &  .or.(lec.and.(.not.loc_lec)) .or. (lnc.and.(.not.loc_lnc))  & 
     &    ) then
          if(mype.eq.0) then
            write(*,*) 'Paramesh Error: for a call to', & 
     &       '  mpi_amr_1blk_guardcell', & 
     &       ' one of more of the arguments lcc/lfc/lec/lnc are not', & 
     &       ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
            print *,' nvar = ',nvar,lcc,loc_lcc
          endif
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
        endif

      else

        loc_lcc = .false.
        if(nvar_work.gt.0)     loc_lcc = .true.
        if(lcc.and.(.not.loc_lcc)) then
          if(mype.eq.0) then
            write(*,*) 'Paramesh Error: for a call to', & 
     &       '  mpi_amr_1blk_guardcell', & 
     &       ' one of more of the arguments lcc/lfc/lec/lnc are not', & 
     &       ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
            print *,' nvar = ',nvar,lcc,loc_lcc
          endif
          call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
        endif

      end if

      end if

#ifdef DEBUG
      write(*,*) 'amr_1blk_guardcell: pe ',mype,' err check ended'
#endif /* DEBUG */

      surrblks  = -1          ! catchy initialization
      tsurrblks = -1  

! construct a list of blocks surrounding local block lb
      l_parent = .true.
      if(l_srl_only) l_parent = .false.
      call mpi_amr_local_surr_blks_lkup(mype,lb, & 
     &                          surrblks,l_parent,psurrblks)

! relate surrblks with the guard block indices stored implicitly 
! in laddress and update tsurrblks

      tsurrblks = surrblks    ! set up temporary array
#ifdef DEBUG
       if(mype.eq.0) & 
     &        write(*,*) '1blk_guardcell : pe ',mype, & 
     &         ' working on lb ',lb,' tsurrblks ', & 
     &         tsurrblks(:,:,:,2), & 
     &         ' lfc ',lfc
#endif /*DEBUG */


      if (timing_mpi) then
         time2 = mpi_wtime()
      endif
! guard block indeces start at strt_buffer after lnblocks, and end at
! last_buffer as determined in subroutine mpi_commatrix.


       do k = 2-k3d,2+k3d     ! loop over all its surrounding blocks
       do j = 2-k2d,2+k2d
       do i = 1,3
         if( surrblks(2,i,j,k).ne.mype .and. & 
     &       surrblks(2,i,j,k).ge.0 ) then

!-
         lfound = .false.
         ll = ladd_strt(surrblks(2,i,j,k))
         do while(.not.lfound.and.ll.le.ladd_end(surrblks(2,i,j,k)))
!         ll = strt_buffer
!         do while(.not.lfound.and.ll.le.last_buffer)
           if( (surrblks(2,i,j,k).eq.laddress(2,ll))  .and. & 
     &         (surrblks(1,i,j,k).eq.laddress(1,ll)) ) then
                              ! found the corresponding block id ll
             tsurrblks(1,i,j,k) = ll
             tsurrblks(2,i,j,k) = mype
             tsurrblks(3,i,j,k) = nodetype(ll)     !?????
             lfound = .true.
#ifdef DEBUG
       if(mype.eq.0) & 
     &       write(*,*) 'pe ',mype,' looking for ',surrblks(:,i,j,k) & 
     &             ,' in slot ',ll,' FOUND ', & 
     &             laddress(:,ll)
#endif /*DEBUG */
           else
#ifdef DEBUG
       if(mype.eq.0) & 
     &       write(*,*) 'pe ',mype,' looking for ',surrblks(:,i,j,k) & 
     &             ,' in slot ',ll,' found ', & 
     &             laddress(:,ll)
#endif /*DEBUG */
             ll = ll+1  
           endif
         enddo
!-

         endif
         if( (tsurrblks(2,i,j,k).ne.mype) .and. & 
     &       (tsurrblks(2,i,j,k).ne.-1) .and. & 
     &       (tsurrblks(1,i,j,k).gt.-20) ) then
             write(*,*) 'ERROR in mpi_amr_1blk_guardcell : pe ',mype, & 
     &         ' working on lb ',lb,' neigh ',i,j,k, & 
     &         ' cannot find surrblk ', & 
     &         surrblks(:,i,j,k),' on this proc ', & 
     &         ' laddress ',laddress(:,strt_buffer:last_buffer), & 
     &         ' strt_buffer,last_buffer ',strt_buffer,last_buffer, & 
     &         ' tsurrblks ',tsurrblks(1:2,i,j,k) & 
     &       ,' ladd_strt ',ladd_strt,' ladd_end ',ladd_end & 
     &       ,' comm pattern id ',mpi_pattern_id
             call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
         endif
        enddo
        enddo
        enddo


! update surrblks with the local guard block info

      surrblks = tsurrblks

#ifdef DEBUG
       if(mype.eq.0) then
        print *, 'MPI SURRBLKS on proc=',mype,' and block=',lb
        do j = 2+k2d,2-k2d,-1
        write(*,'(i4,2x,i4,2x,i4)') (surrblks(1,i,j,2),i=1,3) 
        enddo
       endif
#endif

       if(spherical_pm) then
       ipolar = 0
       if(lsingular_line) then
       if(abs(bnd_box(1,2,lb)).lt.eps) ipolar(1) = -1
       if(abs(bnd_box(2,2,lb)-pi).lt.eps) ipolar(2) = 1
       endif
       endif

      if (timing_mpi) then
      timer_amr_1blk_guardcell(1) = timer_amr_1blk_guardcell(1) & 
     &                           +  mpi_wtime() - time2
      time2 = mpi_wtime()
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(.not.l_srl_only) then
!
! Are there any coarse neighbors?
        lcoarse = .false.


        if(parent(1,lb).gt.0) then
        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
          if(surrblks(1,i,j,k).gt.-20.and.surrblks(1,i,j,k).lt.0) & 
     &                   lcoarse = .true.
        enddo
        enddo
        enddo
        endif 

        endif                        ! end of l_srl_only if test

!------------------------------------
!
! Put leaf block lb's data into the data_1blk.fh datastructures, 
! with the appropriate guardcell padding.

          idest = 1

          call amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)

#ifdef DEBUGY
       if(mype.eq.0) then
          write(*,*) 'after perm : pe ',mype,' blk ',lb,pe, & 
     &    ' facevarx1 ',facevarx1(1,:,:,1,1)
       endif
#endif /* DEBUG */

!------------------------------------
!
        if(iopt.eq.1) then
          pcache_pe  = pcache_pe_u
          pcache_blk = pcache_blk_u
        elseif(iopt.ge.2) then
          pcache_pe  = pcache_pe_w
          pcache_blk = pcache_blk_w
        endif

#ifdef DEBUG
        write(*,*) 'pe ',mype,'blk ',lb, & 
     &      'lcoarse ',lcoarse,' l_srl_only ',l_srl_only
#endif

        if(.not.l_srl_only) then
        if(lcoarse) then

!
! Put data from lb's parent into the data_1blk.fh datastructures, with the
! appropriate guardcell padding. Check to see if data is currently cached.
          parent_lb = parent(1,lb)
          parent_pe = parent(2,lb)


#ifdef DEBUG
       if(mype.eq.0) then
        write(*,*) 'pe ',mype,'blk ',lb, & 
     &      'parent blk ',parent_lb,parent_pe,' cached ', & 
     &      pcache_blk,pcache_pe
        endif
#endif
        if( (parent_lb.gt.0) .and. & 
     &      ((parent_lb.ne.pcache_blk).or.(parent_pe.ne.pcache_pe) ) & 
     &      ) then

! record id of new parent block placed in cache
          lnew_parent = .true.
          pcache_blk = parent_lb
          pcache_pe  = parent_pe

#ifdef DEBUG
       if(mype.eq.0) then
        write(*,*) 'pe ',mype,'blk ',lb, & 
     &      'fetching parent blk ',parent_lb,parent_pe, & 
     &      'new cache ',pcache_blk,pcache_pe,' lnew_parent ', & 
     &      lnew_parent,' parent(:,lb) ',parent(:,lb)
        endif
#endif

        if(lcc) then
           if (first_cc) then
              unk1(:,:,:,:,2) = 0.
              work1(:,:,:,2) = 0.
              first_cc = .false.
           end if
        end if
        if(lnc) then
           if (first_nc) then
              unk_n1(:,:,:,:,2) = 0.
              first_nc = .false.
           end if
        end if
        if(lec) then
           if (first_ec) then
             unk_e_x1(:,:,:,:,2) = 0.
             unk_e_y1(:,:,:,:,2) = 0.
             unk_e_z1(:,:,:,:,2) = 0.
             first_ec = .false.
           endif
        endif
        if(lfc) then
           if (first_fc) then
             facevarx1(:,:,:,:,2) = 0.
             facevary1(:,:,:,:,2) = 0.
             facevarz1(:,:,:,:,2) = 0.
             first_fc = .false.
           endif
        end if

#ifdef DEBUG
        write(*,*) 'pe ',mype,'blk ',lb,'calling getremoteblock'
#endif
        idest = 2

        call mpi_amr_get_remote_block(mype,parent_pe,parent_lb, & 
     &                                idest,iopt,lcc,lfc,lec,lnc, & 
     &                            nlayers0x,nlayers0y,nlayers0z)

!
!------------------------------------
! Do guardcell filling for lb's parent from any surrounding blocks at 
! the same refinement level as this parent.
! Diagonal elements are required to ensure that all cells are filled
! correctly when icoord is non-zero.
          iblock=2
          icoord_loc = 0
          ldiag_loc = .true.

       if(spherical_pm) then
       if(parent_pe.eq.mype) then
         pbnd_box(:,2) = bnd_box(:,2,parent_lb)
       else
         ijk = 1
         do ij=strt_buffer,last_buffer
           if(laddress(1,ij).eq.parent_lb & 
     &     .and.laddress(2,ij).eq.parent_pe) then
              ijk = ij
              pbnd_box(:,2) = bnd_box(:,2,ijk)
           endif
         enddo
         if(ijk.eq.-1) then
           write(*,*)  & 
     &      'mpi_amr_1blk_guardcell : ijk still -1 : search failed'
           call amr_abort()
         endif
       endif
       ippolar = 0
       if(lsingular_line) then
       if(abs(pbnd_box(1,2)).lt.eps) ippolar(1) = -1
       if(abs(pbnd_box(2,2)-pi).lt.eps) ippolar(2) = 1
       endif
       endif

#ifdef DEBUG
        write(*,*) 'pe ',mype,'blk ',lb,'calling gsrl'
#endif
          call amr_1blk_guardcell_srl(mype,parent_pe,parent_lb, & 
     &                                iblock,iopt,nlayers,psurrblks, & 
     &                                lcc,lfc,lec,lnc, & 
     &                                icoord_loc,ldiag_loc, & 
     &                     nlayers0x,nlayers0y,nlayers0z,ippolar)

         if (lcc .and. iopt.eq.1)  & 
     &        call flash_convert_cc_hook(unk1(:,:,:,:,2), nvar, & 
     &         il_bnd1,iu_bnd1, jl_bnd1,ju_bnd1, kl_bnd1,ku_bnd1, & 
     &         why=gr_callReason_PROLONG)


        endif       ! end if parents data not previously cached

!
!------------------------------------
!
! Do guardcell filling from coarse neigbors into the current block
#ifdef DEBUG
        write(*,*) 'pe ',mype,'blk ',lb,'calling gc2f'
#endif
        call mpi_amr_1blk_guardcell_c_to_f( mype,lb,pe,iopt,nlayers, & 
     &                                      surrblks, & 
     &                                      lcc,lfc,lec,lnc, & 
     &                                      icoord,ldiag, & 
     &                                      nlayers0x, & 
     &                                      nlayers0y, & 
     &                                      nlayers0z,ipolar)

        if (lcc .and. iopt.eq.1)  & 
     &       call flash_unconvert_cc_hook(unk1(:,:,:,:,1), nvar, & 
     &       il_bnd1,iu_bnd1, jl_bnd1,ju_bnd1, kl_bnd1,ku_bnd1, & 
     &       where=gr_cells_GUARD, why=gr_callReason_PROLONG)
!------------------------------------

        endif                       ! end of lcoarse if test

        endif                       ! end of l_srl_only if test

      if (timing_mpi) then
      timer_amr_1blk_guardcell(2) = timer_amr_1blk_guardcell(2) & 
     &                           +  mpi_wtime() - time2
      time2 = mpi_wtime()
      endif
!
!------------------------------------
!
! Do guardcell filling from any surrounding blocks at the same refinement
! level as block lb.
        iblock = 1
#ifdef DEBUG
          if(iopt.gt.1)  & 
     &     write(*,*) 'in guardcell_srl : pe ',mype,' blk ',lb, & 
     &        ' surrblks(:,:,:,2) ',surrblks(:,:,:,2)
          if(iopt.eq.1.and.lb.eq.1)  write(*,*) 'located srl call '
#endif /*DEBUG */

        call amr_1blk_guardcell_srl(mype,mype,lb, & 
     &                              iblock,iopt,nlayers,surrblks, & 
     &                              lcc,lfc,lec,lnc, & 
     &                              icoord,ldiag, & 
     &                              nlayers0x,nlayers0y,nlayers0z, & 
     &                              ipolar)


      if (timing_mpi) then
      timer_amr_1blk_guardcell(3) = timer_amr_1blk_guardcell(3) & 
     &                           +  mpi_wtime() - time2
      endif

      if (.not.advance_all_levels) then
!        call amr_1blk_guardcell_f_to_c(mype,pe,lb,
!     .                              iblock,iopt,nlayers,surrblks,
!     .                              lcc,lfc,lec,lnc,icoord,ldiag,
!     .                              nlayers0x,nlayers0y,nlayers0z)
      end if

#ifdef DEBUG
          write(*,*) 'in guardcell_srl : pe ',mype,' blk ',lb, & 
     &        ' surrblks(:,:,:,2) ',surrblks(:,:,:,2)
#endif /*DEBUG */
    
#ifdef DEBUGY
          write(*,*) 'after guardcell_srl : pe ',mype,' blk ',lb, & 
     &    ' unk1 ',unk1(1,:,:,1,1)
#endif /*DEBUG */

        if(iopt.eq.1) then
          pcache_pe_u  = pcache_pe
          pcache_blk_u = pcache_blk
        elseif(iopt.ge.2) then
          pcache_pe_w  = pcache_pe
          pcache_blk_w = pcache_blk
        endif

      if (timing_mpi) then
      timer_amr_1blk_guardcell(0) = timer_amr_1blk_guardcell(0) & 
     &                           +  mpi_wtime() - time1
      endif

#ifdef DEBUG_FLOW_TRACE
      write(*,*)'exiting amr_1blk_guardcell: pe ',mype,' exiting', & 
     &          ' iopt ',iopt
#endif /* DEBUG_FLOW_TRACE */

      return
      end subroutine amr_1blk_guardcell



