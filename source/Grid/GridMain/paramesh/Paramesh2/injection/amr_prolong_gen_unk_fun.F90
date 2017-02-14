!!****if* source/Grid/GridMain/paramesh/Paramesh2/injection/amr_prolong_gen_unk_fun
!!
!! NAME
!!
!!  amr_prolong_gen_unk_fun
!!
!!
!! SYNOPSIS
!!
!!  amr_prolong_gen_unk_fun(recv, ia, ib, ja, jb, ka, kb, isg, &
!!                          ioff, joff, koff, mype)
!!
!!
!!  amr_prolong_gen_unk_fun(real(), integer, integer, integer, integer, &
!!                          integer, integer, integer, &
!!                          integer, integer, integer, integer)
!!
!!
!! DESCRIPTION
!!
!!
!!  This routine takes data from the array recv, originally extracted 
!!  from the solution array unk, and performs a prolongation operation 
!!  on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!  The data in recv is from a parent block (and may have come from off
!!  processor) and the result of the prolongation operation is written 
!!  directly into block isg, which is one of its children (and on the 
!!  local processor).  The position of the child within the parent block 
!!  is specified by the ioff, joff, and koff arguments (given in terms 
!!  of parent block zones, excluding guardcells).
!!
!!  This particular prolongation is conservative, cell-averaged, direct
!!  copy.  This is first order accurate only.  Because it is direct copy,
!!  it is accurate for any geometry.
!!
!!
!! ARGUMENTS
!!
!!  recv          the coarse parent data that will be used to initialize the
!!                child block.  An entire block (including guardcells) of 
!!                coarse data is provided here
!!
!!  ia, ib        the range of *child* zones to fill via prolongation. 
!!  ja, jb        (sometimes we are only doing the guardcells).  This is set 
!!  ka, kb        by amr_prolong_unk_fun, the wrapper for the generic 
!!                prolongation.
!!
!!  isg           the block number of the child block whose data we are 
!!                filling.
!!
!!  ioff          the offsets of the child block into the parent (in terms of 
!!  joff          parent zones, excluding guardcells).
!!  koff
!!
!!  mype          the local processor number (why is this here?)
!!
!!
!! NOTES
!!
!!  this is direct copy interpolantion.
!!
!!  no monotonicity hacks are required, as linear interpolation guarantees that
!!  the prolonged data falls between the coarse cell values.
!!
!!***

subroutine amr_prolong_gen_unk_fun(recv, ia, ib, ja, jb, ka, kb, isg, & 
     ioff, joff, koff, mype, isrc)

  use dBase, ONLY : nxb, nyb, nzb, k2d, k3d, nguard, nvar, ndim, &
       iLo_gc, iHi_gc, jLo_gc, jHi_gc, kLo_gc, kHi_gc, ionmax, &
       geom_cartesian, geom_cylrad, geom_sphrad, &
       geom_cylang, geom_sphtheta, geom_sphphi, &
       CARTESIAN, CYLINDRICAL, SPHERICAL, &
       dBaseGetDataPtrSingleBlock,     &
       GC,                             &
       dBaseReleaseDataPtrSingleBlock, &
       dBaseGetCoords, dBaseKeyNumber, dBaseSpecies

  use runtime_parameters

  implicit none


  real :: recv(nvar,iLo_gc:iHi_gc,jLo_gc:jHi_gc,kLo_gc:kHi_gc)

  integer :: ia, ib, ja, jb, ka, kb

  integer :: isg
  integer :: ioff, joff, koff

  integer :: mype
  integer,intent(IN),OPTIONAL :: isrc

  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1

  integer :: ivar
  integer :: i, j, k, n

  integer :: ico, jco, kco

  integer, save :: inuc_begin

  real :: error

  real, DIMENSION(:,:,:,:), POINTER :: unk

  logical, save :: firstCall = .true.


  if (firstCall) then

     inuc_begin = dBaseSpecies(1)

     firstCall = .false.

  endif

  unk => dBaseGetDataPtrSingleBlock(isg, GC)

! loop over zones in the target (fine/child) block, and prolong the coarse
! data to each child zone.  
!
! (i,j,k) are the indices of the fine zones.  
!
! (ip,jp,kp) are the indices of the coarse zone enclosing each fine zone.  
!
! ip1 and im1 refer to coarse zones offset by one to the right and left, 
! respectively; likewise for jp1, jm1, kp1, km1.  
!
! ico, jco, and kco indicate the location of the child in the parent zone in
! each coordinate direction.

  do k = ka, kb

     kp  = ((k-1)/2) + 1 + (nguard/2)*k3d + koff
     kp1 = kp + k3d
     km1 = kp - k3d

! kco tells us where the child falls in the coarse zone.  kco = 1 means it 
! shares the lower z interface, kco = 0 means it shares the upper z interface 
! of the parent zone
     kco = mod(k, 2)

     do j = ja, jb

        jp  = ((j-1)/2) + 1 + (nguard/2)*k2d + joff
        jp1 = jp + k2d
        jm1 = jp - k2d

        jco = mod(j, 2)

        do i = ia, ib

           ip  = ((i-1)/2) + 1 + (nguard/2) + ioff
           ip1 = ip + 1
           im1 = ip - 1

           ico = mod(i, 2)

           do ivar = 1, nvar

              unk(ivar,i,j,k) = recv(ivar,ip,jp,kp)

           enddo
            
        enddo
     enddo
  enddo

  call dBaseReleaseDataPtrSingleBlock(isg, unk)

  return
end subroutine amr_prolong_gen_unk_fun

