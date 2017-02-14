

!******************************************************************************

!  Routine:     ref_marking()

!  Description: Compute the refinement criterion on all blocks and use it to
!               determine whether each block should be refined or derefined.

!  Written:     K. Olson 8/97, based on Lohner (1987), Comput. Meth. in Appl.
!               Mech. and Eng., vol. 61, pp. 323-338, and on notes and code
!               written by R. deFainchtein.

!  Inputs:      ctore           - Refinement cutoff.  Refine only if
!                                 estimator is greater than this value.
!                                 Should be between 0 and 1.
!               ctode           - Derefinement cutoff.  Derefine only if
!                                 estimator is smaller than this value.
!                                 Should be between 0 and 1.
!               epsil           - Filter coefficient used in computing the
!                                 second-derivative estimator.  See the
!                                 User's Guide for details.
!               lref_max        - Maximum allowed level of refinement.
!               iref            - Index into unk() of variable to refine on.

!  Outputs:     Blocks are marked appropriately (using the global variables
!               refine() and derefine()).  Refinements and derefinements are
!               carried out by amr_refine_derefine().


subroutine ref_marking (ctore, ctode, epsil, lref_max, msgbuffer, iref)

!==============================================================================

  use physicaldata
  use tree
  use Grid_data, ONLY: gr_geometry, gr_oneBlock, gr_meshComm
  implicit none
  include 'mpif.h'

#include "constants.h"

  logical msgbuffer
  integer ndim2
  parameter (ndim2=ndim*ndim)

  logical first

  real delx,dely,delz
  real dely_f, delz_f
  real delu(ndim,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
  real delua(ndim,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

  real delu2(ndim2), delu3(ndim2), delu4(ndim2)

  real num,denom,epsil,error(maxblocks),error_par(maxblocks), & 
       &     errort
  real ctore,ctode,J_sq,xx,dd
  real error_max
  real vx,vy,vz,bv,bx,by,bz,bx1,by1,bz1,rhbi,rhi,bk,et
  real px,py,pz,rho,b2,v2,gf,vv,gr,dx,vt,press
  real times,time_exe

  integer lb,i,j,k,kk,lref_max,iref,nref
  integer mype
  integer ineigh,ineigh_proc,ierr
  integer kstart,kend,jstart,jend,istart,iend
  integer nsend,nrecv
  integer reqr(maxblocks),reqs(maxblocks_tr)
  integer statr(MPI_STATUS_SIZE,maxblocks)
  integer stats(MPI_STATUS_SIZE,maxblocks_tr)
!
!       for the message buffering
!
  integer stags(maxblocks_tr)
  real    svals(maxblocks_tr)
  integer sprocs(maxblocks_tr)
  integer rtags(maxblocks_tr)
  real    rvals(maxblocks_tr)
  integer rprocs(maxblocks_tr)

  real, parameter  :: epsilon = TINY(1.0)

!==============================================================================


  call MPI_Comm_Rank (gr_meshComm, mype, ierr)

#define XCOORD(I) (gr_oneBlock(lb)%firstAxisCoords(CENTER,I))
#define YCOORD(I) (gr_oneBlock(lb)%secondAxisCoords(CENTER,I))

  do lb = 1,lnblocks

     error(lb) = 0.

     if (nodetype(lb).eq.1.or.nodetype(lb).eq.2) then

#if N_DIM >= 2
!            dely = bsize(2,lb)/float(nyb)
!            dely = 1./(2.*dely)

        dely = 0.5e0*float(nyb)/bsize(2,lb)
        dely_f = dely
#endif

#if N_DIM == 3
!            delz = bsize(3,lb)/float(nzb)
!            delz = 1./(2.*delz)

        delz = 0.5e0*float(nzb)/bsize(3,lb)
        delz_f = delz
#endif

! Compute first derivatives

        do k = 1+k3d,nzb+(k3d*((2*nguard)-1))
           do j = 1+k2d,nyb+(k2d*((2*nguard)-1))
              do i = 2,nxb+(2*nguard)-1

                 delx = 1.0/(XCOORD(i+1) - XCOORD(i-1))
                 
                 ! d/dx                
                 delu(1,i,j,k) = unk(iref,i+1,j,k,lb) - unk(iref,i-1,j,k,lb)
                 delu(1,i,j,k) = delu(1,i,j,k)*delx
                 
                 delua(1,i,j,k) = abs(unk(iref,i+1,j,k,lb)) + &
                      abs(unk(iref,i-1,j,k,lb))
                 delua(1,i,j,k) = delua(1,i,j,k)*delx

#if N_DIM >= 2
                 ! d/dy
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    dely_f = dely/XCOORD(i)
                 end if

                 delu(2,i,j,k) = unk(iref,i,j+1,k,lb) - unk(iref,i,j-1,k,lb)
                 delu(2,i,j,k) = delu(2,i,j,k)*dely_f

                 delua(2,i,j,k) = abs(unk(iref,i,j+1,k,lb)) +  & 
                      abs(unk(iref,i,j-1,k,lb))
                 delua(2,i,j,k) = delua(2,i,j,k)*dely_f
#endif


#if N_DIM == 3
                 ! d/dz
                 if (gr_geometry == SPHERICAL) then
                    delz_f = delz/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    delz_f = delz/XCOORD(i)
                 end if
                 delu(3,i,j,k) = unk(iref,i,j,k+1,lb) - unk(iref,i,j,k-1,lb)
                 delu(3,i,j,k) = delu(3,i,j,k)*delz_f

                 delua(3,i,j,k) = abs(unk(iref,i,j,k+1,lb)) +  & 
                      abs(unk(iref,i,j,k-1,lb))
                 delua(3,i,j,k) = delua(3,i,j,k)*delz_f
#endif

              end do
           end do
        end do


! Compute second derivatives

! Test if at a block boundary

! Two guardcells            
        kstart = 2*k3d+1
        kend   = nzb+(k3d*((2*nguard)-2))

        jstart = 2*k2d+1
        jend   = nyb+(k2d*((2*nguard)-2))

        istart = 3
        iend   = nxb+(2*nguard)-2

! One guardcell
!            kstart = 2*k3d+1+k3d
!            kend   = nzb+(k3d*((2*nguard)-2))-k3d
!            jstart = 2*k2d+1+k2d
!            jend   = nyb+(k2d*((2*nguard)-2))-k2d
!            istart = nguard
!            iend   = nxb+(2*nguard)-3
! No guardcells
!            kstart = k3d*nguard+1
!            kend   = nzb+k3d*nguard
!            jstart = k2d*nguard+1
!            jend   = nyb+k2d*nguard
!            istart = nguard+1
!            iend   = nxb+nguard

! only compute the derivatives inside the physical domain -- don't venture
! into the boundary conditions

        if (neigh(1,1,lb).le.-20) istart = nguard+1
        if (neigh(1,2,lb).le.-20) iend   = nguard+nxb

#if N_DIM >= 2
        if (neigh(1,3,lb).le.-20) jstart = nguard*k2d+1
        if (neigh(1,4,lb).le.-20) jend   = nguard*k2d+nyb
#endif

#if N_DIM == 3
        if (neigh(1,5,lb).le.-20) kstart = nguard*k3d+1
        if (neigh(1,6,lb).le.-20) kend   = nguard*k3d+nzb
#endif

        do k = kstart,kend
           do j = jstart,jend
              do i = istart,iend


                 delx = 1.0/(XCOORD(i+1) - XCOORD(i-1))

! d/dxdx
                 delu2(1) = delu(1,i+1,j,k) - delu(1,i-1,j,k)
                 delu2(1) = delu2(1)*delx

                 delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
                 delu3(1) = delu3(1)*delx

                 delu4(1) = delua(1,i+1,j,k) + delua(1,i-1,j,k)
                 delu4(1) = delu4(1)*delx

#if N_DIM >= 2
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    dely_f = dely/XCOORD(i)
                 end if

! d/dydx
                 delu2(2) = delu(1,i,j+1,k) - delu(1,i,j-1,k)
                 delu2(2) = delu2(2)*dely_f

                 delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
                 delu3(2) = delu3(2)*dely_f

                 delu4(2) = delua(1,i,j+1,k) + delua(1,i,j-1,k)
                 delu4(2) = delu4(2)*dely_f

! d/dxdy
                 delu2(3) = delu(2,i+1,j,k) - delu(2,i-1,j,k)
                 delu2(3) = delu2(3)*delx

                 delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
                 delu3(3) = delu3(3)*delx

                 delu4(3) = delua(2,i+1,j,k) + delua(2,i-1,j,k)
                 delu4(3) = delu4(3)*delx

! d/dydy
                 delu2(4) = delu(2,i,j+1,k) - delu(2,i,j-1,k)
                 delu2(4) = delu2(4)*dely_f

                 delu3(4) = abs(delu(2,i,j+1,k)) + abs(delu(2,i,j-1,k))
                 delu3(4) = delu3(4)*dely_f

                 delu4(4) = delua(2,i,j+1,k) + delua(2,i,j-1,k)
                 delu4(4) = delu4(4)*dely_f
#endif

#if N_DIM == 3
                 if (gr_geometry == SPHERICAL) then
                    delz_f = delz/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    delz_f = delz/XCOORD(i)
                 end if

! d/dzdx
                 delu2(5) = delu(1,i,j,k+1) - delu(1,i,j,k-1)
                 delu2(5) = delu2(5)*delz_f

                 delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
                 delu3(5) = delu3(5)*delz_f

                 delu4(5) = delua(1,i,j,k+1) + delua(1,i,j,k-1)
                 delu4(5) = delu4(5)*delz_f

! d/dzdy
                 delu2(6) = delu(2,i,j,k+1) - delu(2,i,j,k-1)
                 delu2(6) = delu2(6)*delz_f

                 delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
                 delu3(6) = delu3(6)*delz_f

                 delu4(6) = delua(2,i,j,k+1) + delua(2,i,j,k-1)
                 delu4(6) = delu4(6)*delz_f

! d/dxdz
                 delu2(7) = delu(3,i+1,j,k) - delu(3,i-1,j,k)
                 delu2(7) = delu2(7)*delx

                 delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
                 delu3(7) = delu3(7)*delx

                 delu4(7) = delua(3,i+1,j,k) + delua(3,i-1,j,k)
                 delu4(7) = delu4(7)*delx

! d/dydz
                 delu2(8) = delu(3,i,j+1,k) - delu(3,i,j-1,k)
                 delu2(8) = delu2(8)*dely_f

                 delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
                 delu3(8) = delu3(8)*dely_f

                 delu4(8) = delua(3,i,j+1,k) + delua(3,i,j-1,k)
                 delu4(8) = delu4(8)*dely_f

! d/dzdz
                 delu2(9) = delu(3,i,j,k+1) - delu(3,i,j,k-1)
                 delu2(9) = delu2(9)*delz_f

                 delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
                 delu3(9) = delu3(9)*delz_f

                 delu4(9) = delua(3,i,j,k+1) + delua(3,i,j,k-1)
                 delu4(9) = delu4(9)*delz_f
#endif

! compute the error
                 num = 0.
                 denom = 0.

                 do kk = 1, ndim2
                    num = num + delu2(kk)**2
                    denom = denom + (delu3(kk) + (epsil*delu4(kk)+epsilon))**2
                 end do

! mz -- compare the square of the error 
                 error(lb) = max(error(lb),num/denom)

              end do
           end do
        end do

! store the maximum error for the current block
        error(lb) = sqrt(error(lb))

     end if

  end do

! MARK FOR REFINEMENT OR DEREFINEMENT

! first communicate error of parent to its children
! parents collect messages from children

  error_par(1:lnblocks) = 0.
  nrecv = 0
  do lb = 1,lnblocks
     if(parent(1,lb).gt.-1) then
        if (parent(2,lb).ne.mype) then
           nrecv = nrecv + 1
           if (msgbuffer) then
              rprocs(nrecv) = parent(2,lb)
              rtags (nrecv) = lb
           else
              call MPI_IRecv(error_par(lb),    & 
                   &                        1, & 
                   &                        MPI_DOUBLE_PRECISION, & 
                   &                        parent(2,lb), & 
                   &                        lb, & 
                   &                        gr_meshComm, & 
                   &                        reqr(nrecv), & 
                   &                        ierr)
           endif
        else
           error_par(lb) = error(parent(1,lb))
        end if
     end if
  end do

! parents send error to children

  nsend = 0
  do lb = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,lb).gt.-1) then
           if (child(2,j,lb).ne.mype) then
              nsend = nsend + 1
              if (msgbuffer) then
                 sprocs(nsend) = child(2,j,lb)
                 stags (nsend) = child(1,j,lb)
                 svals (nsend) = error(lb)
              else 
                 call MPI_ISend(error(lb), & 
                      &                           1, & 
                      &                           MPI_DOUBLE_PRECISION, & 
                      &                           child(2,j,lb), &  ! PE TO SEND TO
                      &                           child(1,j,lb), &  ! THIS IS THE TAG
                      &                           gr_meshComm, & 
                      &                           reqs(nsend), & 
                      &                           ierr)
              endif
           end if
        end if
     end do
  end do

  if (.not.(msgbuffer)) then
     if (nsend.gt.0) then
        call MPI_Waitall (nsend, reqs, stats, ierr)
     end if
     if (nrecv.gt.0) then
        call MPI_Waitall (nrecv, reqr, statr, ierr)
     end if
  else 
     call b_dbl_sendrcv(1020,1,nsend,sprocs,stags,svals, & 
          &                            nrecv,rprocs,rtags,rvals)
     nrecv = 0
     do lb = 1,lnblocks
        if(parent(1,lb).gt.-1) then
           if (parent(2,lb).ne.mype) then
              nrecv = nrecv + 1
              error_par(lb) = rvals(nrecv)
           endif
        endif
     enddo
  endif


  do lb = 1,lnblocks

     if (nodetype(lb).eq.1) then


! test for derefinement

        if (.not.refine(lb).and..not.stay(lb) & 
             &          .and.error(lb).le.ctode & 
             &          .and.error_par(lb).le.ctode) then
           derefine(lb) = .TRUE.
        else
           derefine(lb) = .FALSE.
        end if

! test for refinement

        if (error(lb).gt.ctore) then
           derefine(lb) = .FALSE.
           refine(lb) = .TRUE.
        end if

        if (error(lb).gt.ctode.or.error_par(lb).gt.ctode)  & 
             &           stay(lb) = .TRUE.

        if (lrefine(lb).ge.lref_max)  & 
             &           refine(lb) = .FALSE.

!        print *, lb, error(lb), error_par(lb), &
!     &          refine(lb), derefine(lb), stay(lb)

     end if

  end do

!==============================================================================

  return
end subroutine ref_marking


