       subroutine amr_guardcell_cc_srl(mype,iopt,nlayers,idiag,idir,maxNodetype_gcWanted)


! $RCSfile: amr_guardcell_cc_srl.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks assuming that exchange is only required between blocks at the
! same refinement level. This routine operates on data defined at cell
! centers (ie. unk or work).
! The actual exchanges are performed with calls
! to the routines cp_loc and cp_remote.
!
! The order of the loops (ie placing the loop over faces on the outside,
! with a barrier after each iteration) is required to guarantee 
! that the guardcells along block edges are correct.
!
!
! Written :     Peter MacNeice          December 1996
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor number
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt=2 will use 'work'
!      nlayers        the number of guard cell layers at each boundary
!
!------------------------------------
use physicaldata
      use tree
      use workspace
      implicit none
      include         'mpif.h'



      integer mype,iopt,nlayers,idiag,idir
      integer,intent(IN) :: maxNodetype_gcWanted
      integer nodetypet
      integer reqr(maxblocks)
      integer statr(MPI_STATUS_SIZE,maxblocks)
      integer block_point(maxblocks)
      integer block_int3d,block_int1,size_of_double
      integer block_int2d,myblockint

      integer j

      integer js,je,jsend,jneigh
      integer istart,iend,jstart,jend,kstart,kend
      integer neighr,neighs,isg
      integer ierr

      real time_comm,time_rest

      save nodetypet

!------------------------------------

      time_comm = 0.
      time_rest = 0.

! All guardcells need to be exchanged here to get correct values in c_to_f
      js = 1
      je = nfaces

!      if (idir.eq.1) then
!        js = 1
!        je = 2
!      elseif (idir.eq.2) then
!        js = 3
!        je = 4
!      elseif (idir.eq.3) then
!        js = 5
!        je = 6
!      end if
!      if (idiag.eq.1) then
!        js = 1
!        je = nfaces
!      end if

! cycle through block faces

      do j = js,je

        if (j.eq.1) jsend = 2
        if (j.eq.2) jsend = 1
        if (j.eq.3) jsend = 4
        if (j.eq.4) jsend = 3
        if (j.eq.5) jsend = 6
        if (j.eq.6) jsend = 5

! cycle through the grid blocks on this processor
        if(lnblocks.gt.0) then

        if (iopt.eq.1) then

        if (j.eq.1.or.j.eq.2) then
          call MPI_TYPE_VECTOR (nyb, & 
     &         nvar*nguard, & 
     &         nvar*iu_bnd, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int2d, & 
     &         ierr)
          myblockint = block_int2d
#if N_DIM == 3
          call MPI_TYPE_HVECTOR (nzb, & 
     &         1, & 
     &         nvar*iu_bnd*ju_bnd*8, & 
     &         block_int2d, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
#endif
          if (j.eq.1) then
            istart = nguard+(nxb-nguard)+1
            jstart = nguard*k2d+1
            kstart = nguard*k3d+1
            iend   = 1
            jend   = nguard*k2d+1
            kend   = nguard*k3d+1
          elseif (j.eq.2) then
            istart = nguard+1
            jstart = nguard*k2d+1
            kstart = nguard*k3d+1
            iend   = nguard+nxb+1
            jend   = nguard*k2d+1
            kend   = nguard*k3d+1
          end if
        end if

        if (j.eq.3.or.j.eq.4) then
          call MPI_TYPE_VECTOR (nguard, & 
     &         nvar*iu_bnd, & 
     &         nvar*iu_bnd, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int2d, & 
     &         ierr)
          myblockint = block_int2d
#if N_DIM == 3
          call MPI_TYPE_HVECTOR (nzb, & 
     &         1, & 
     &         nvar*iu_bnd*ju_bnd*8, & 
     &         block_int2d, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
#endif
          if (j.eq.3) then
            istart = 1
            jstart = (nguard+(nyb-nguard))*k2d+1
            kstart = nguard*k3d+1
            iend   = 1
            jend   = 1
            kend   = nguard*k3d+1
          elseif (j.eq.4) then
            istart = 1
            jstart = nguard*k2d+1
            kstart = nguard*k3d+1
            iend   = 1
            jend   = (nguard+nyb)*k2d+1
            kend   = nguard*k3d+1
          end if
        end if

        if (j.eq.5.or.j.eq.6) then
          call MPI_TYPE_VECTOR(nguard*ju_bnd, & 
     &         nvar*iu_bnd, & 
     &         nvar*iu_bnd, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
          if (j.eq.5) then
            istart = 1
            jstart = 1
            kstart = nguard+(nzb-nguard)+1
            iend   = 1
            jend   = 1
            kend   = 1
          elseif (j.eq.6) then
            istart = 1
            jstart = 1
            kstart = nguard*k3d+1
            iend   = 1
            jend   = 1
            kend   = nguard+nzb+1
          end if
        end if

        if (ndim.eq.3.and.j.ne.5.and.j.ne.6)  & 
     &      call MPI_TYPE_COMMIT(block_int2d,ierr)
        call MPI_TYPE_COMMIT(myblockint,ierr)

        neighr = 0
        do isg = 1,lnblocks
          if(neigh(1,j,isg).gt.-1) then
            if(neigh(2,j,isg).ne.mype) then
              if(nodetype(isg).eq.1.or.nodetype(isg).eq.2 .or. &
                   ((maxNodetype_gcWanted>1) .AND. nodetype(isg).eq.3)) then
                neighr = neighr + 1
! receive guardcells directly into local array
                call MPI_IRECV(unk(1,iend,jend,kend,isg), & 
     &               1, & 
     &               myblockint, & 
     &               neigh(2,j,isg), & 
     &               neigh(1,j,isg), & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
              end if
            end if
          end if
        end do

! send messages if neighbor is off processor

        neighs = 0
        do isg = 1,lnblocks
          if(neigh(1,jsend,isg).gt.-1) then
            if(neigh(2,jsend,isg).ne.mype) then
              if (neigh_type(jsend,isg).eq.1.or. & 
     &            neigh_type(jsend,isg).eq.2 .or. &
     ((maxNodetype_gcWanted>1) .AND. neigh_type(jsend,isg).eq.3)) then
                neighs = neighs + 1
                call MPI_SSEND(unk(1,istart,jstart,kstart,isg), & 
     &               1, & 
     &               myblockint, & 
     &               neigh(2,jsend,isg), &  ! PE TO SEND TO
     &               isg,       &  ! THIS IS THE TAG
     &               MPI_COMM_WORLD, & 
     &               ierr)
              end if
            end if
          end if
        end do

        if (neighr.gt.0) then
          call MPI_WAITALL (neighr, reqr, statr, ierr)
        end if
        
        else                      ! iopt

        if (j.eq.1.or.j.eq.2) then
          call MPI_TYPE_VECTOR (nyb, & 
     &         nguard_work, & 
     &         iuw, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int2d, & 
     &         ierr)
          myblockint = block_int2d
#if N_DIM == 3
          call MPI_TYPE_HVECTOR (nzb, & 
     &         1, & 
     &         iuw*juw*8, & 
     &         block_int2d, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
#endif
          if (j.eq.1) then
            istart = nguard_work+(nxb-nguard_work)+1
            jstart = nguard_work*k2d+1
            kstart = nguard_work*k3d+1
            iend   = 1
            jend   = nguard_work*k2d+1
            kend   = nguard_work*k3d+1
          elseif (j.eq.2) then
            istart = nguard_work+1
            jstart = nguard_work*k2d+1
            kstart = nguard_work*k3d+1
            iend   = nguard_work+nxb+1
            jend   = nguard_work*k2d+1
            kend   = nguard_work*k3d+1
          end if
        end if

        if (j.eq.3.or.j.eq.4) then
          call MPI_TYPE_VECTOR (nguard_work, & 
     &         iuw, & 
     &         iuw, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int2d, & 
     &         ierr)
          myblockint = block_int2d
#if N_DIM == 3
          call MPI_TYPE_HVECTOR (nzb, & 
     &         1, & 
     &         iuw*juw*8, & 
     &         block_int2d, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
#endif
          if (j.eq.3) then
            istart = 1
            jstart = (nguard_work+(nyb-nguard_work))*k2d+1
            kstart = nguard_work*k3d+1
            iend   = 1
            jend   = 1
            kend   = nguard_work*k3d+1
          elseif (j.eq.4) then
            istart = 1
            jstart = nguard_work*k2d+1
            kstart = nguard_work*k3d+1
            iend   = 1
            jend   = (nguard_work+nyb)*k2d+1
            kend   = nguard_work*k3d+1
          end if
        end if

        if (j.eq.5.or.j.eq.6) then
          call MPI_TYPE_VECTOR(nguard_work*juw, & 
     &         iuw, & 
     &         iuw, & 
     &         MPI_DOUBLE_PRECISION, & 
     &         block_int3d, & 
     &         ierr)
          myblockint = block_int3d
          if (j.eq.5) then
            istart = 1
            jstart = 1
            kstart = nguard_work+(nzb-nguard_work)+1
            iend   = 1
            jend   = 1
            kend   = 1
          elseif (j.eq.6) then
            istart = 1
            jstart = 1
            kstart = nguard_work*k3d+1
            iend   = 1
            jend   = 1
            kend   = nguard_work+nzb+1
          end if
        end if

        if (ndim.eq.3.and.j.ne.5.and.j.ne.6)  & 
     &      call MPI_TYPE_COMMIT(block_int2d,ierr)
        call MPI_TYPE_COMMIT(myblockint,ierr)

        neighr = 0
        do isg = 1,lnblocks
          if(neigh(1,j,isg).gt.-1) then
            if(neigh(2,j,isg).ne.mype) then
              if (nodetype(isg).eq.1.or.nodetype(isg).eq.2 .or. &
                   ((maxNodetype_gcWanted>1) .AND. nodetype(isg).eq.3)) then
                neighr = neighr + 1
                call MPI_IRECV(work(iend,jend,kend,isg,1), & 
     &               1, & 
     &               myblockint, & 
     &               neigh(2,j,isg), & 
     &               neigh(1,j,isg), & 
     &               MPI_COMM_WORLD, & 
     &               reqr(neighr), & 
     &               ierr)
              end if
            end if
          end if
        end do

! send messages if neighbor is off processor

        neighs = 0
        do isg = 1,lnblocks
          if(neigh(1,jsend,isg).gt.-1) then
            if(neigh(2,jsend,isg).ne.mype) then
              if (neigh_type(jsend,isg).eq.1.or.&
                  neigh_type(jsend,isg).eq.2.or.&
                  ((maxNodetype_gcWanted>1) .AND. neigh_type(jsend,isg).eq.3)) then
                neighs = neighs + 1
                call MPI_SSEND(work(istart,jstart,kstart,isg,1), & 
     &               1, & 
     &               myblockint, & 
     &               neigh(2,jsend,isg), &  ! PE TO SEND TO
     &               isg,       &  ! THIS IS THE TAG
     &               MPI_COMM_WORLD, & 
     &               ierr)
              end if
            end if
          end if
        end do

        if (neighr.gt.0) then
          call MPI_WAITALL (neighr, reqr, statr, ierr)
        end if

        end if ! iopt
        
 25     do isg = 1,lnblocks

! if current block has no children or if it has at least one leaf
! child then get the guard cell info that it needs.
!!$           if(nodetype(isg).eq.3 .AND. ANY(neigh_type(max(j,jsend)+1:2*ndim,isg)==2)) then
!!$              write(*,*)isg,'j,jsend,neigh_type(j,isg):',j,jsend,neigh_type(j,isg)
!!$              write(*,*)isg,'neigh_type(max(j,jsend)+1:2*ndim,isg):',neigh_type(max(j,jsend)+1:2*ndim,isg)
!!$           end if
        if(nodetype(isg).eq.1.or.nodetype(isg).eq.2 .or. &
             ((maxNodetype_gcWanted>1) .AND. nodetype(isg).eq.3 .AND. ANY(neigh_type(max(j,jsend)+1:2*ndim,isg)==2))) then
          if (neigh(2,j,isg).eq.mype) then
            jneigh = neigh(1,j,isg)
! If neighbor block exists at this refinement level
            if(jneigh.gt.-1) then
              call amr_cp_loc( isg,jneigh,j,iopt,idiag,nlayers )
            end if
          end if
        endif

        enddo

        call MPI_TYPE_FREE(myblockint,ierr)
        if (ndim.eq.3.and.j.ne.5.and.j.ne.6)  & 
     &      call MPI_TYPE_FREE(block_int2d,ierr)

      endif

      enddo

      return
      end
