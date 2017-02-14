!!****if* source/Grid/GridSolvers/Multipole/gr_mpolePotential
!!
!! NAME
!!
!!  gr_mpolePotential
!!
!! SYNOPSIS
!!
!!  gr_mpolePotential(integer, intent(IN): idensvar,
!!                  integer, intent(IN): ipotvar,
!!                  real,intent(IN)    : poisfact)
!!
!! DESCRIPTION
!!
!!  Computes the potential field due to a density field associated
!!  with a given set of multipole moments.  The moments are taken
!!  from the mpole_common variable Moment().  On output, the
!!  variable indexed by ipotvar contains the potential.  This
!!  calculation is entirely local to each processor, as each
!!  processor now has a separate copy of the moments.
!!
!! PARAMETERS
!!
!!  mpole_subSample  -- integer to control number of subzones in each direction
!!                           for smoothing potential calculations.
!!                      Also slows down calculation considerably when != 1
!!
!! NOTES
!!
!!  2-D spherical version adopted after M. Steinmetz, MPA, Garching, Germany.
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePotential (idensvar, ipotvar, poisfact)

  use gr_mpoleData, only : G_1DSPHERICAL,mpole_geometry, fourpi,&
                         Moment,Momtmp,gbnd, lstep,Even, Odd, &
                         qmax, mpole_lmax, mpole_mmax,twopi,&
                         G_2DCYLINDRICAL,G_1DSPHERICAL,G_3DCARTESIAN,G_3DAXISYMMETRIC,&
                         G_2DSPHERICAL,& 
                         Xcm,Ycm,Zcm,mpole_subSample,mpole_subSampleInv, fourpi_inv, &
                         point_mass, Newton
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkBoundBox, Grid_getDeltas, &
    Grid_getBlkPtr, Grid_getBlkIndexLimits, Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshComm

  implicit none

#include "Flash.h"
#include "constants.h"
  include "Flash_mpi.h"

  integer,intent(IN)  :: idensvar, ipotvar
  real,intent(IN)      :: poisfact
  
  integer   :: i, j, k, l,m, lb, error
  real      :: potential, delx, dely, delz, dvol, x, y, z
  real      :: mpfactor
  
  integer   :: ii, jj, kk
  real      :: xx, yy, zz, potsum
  real      :: delxx, delyy, delzz, ddvol

  integer :: blockCount, blockList(MAXBLOCKS)

  real,dimension(MDIM) :: delta
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real, pointer,dimension(:,:,:,:)      :: solnData
  
  integer :: imax, jmax, kmax, imin, jmin, kmin

  real :: gbnd1

  

  !=========================================================================
  
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
     
  if(mpole_geometry /= G_2DSPHERICAL) then

     !Scale the Poisson source term factor appropriately.

     mpfactor = poisfact * fourpi_inv

     !Compute potential on all locally held blocks.


     do lb = 1, blockCount

        !get index size of blk (used to be nxb, nyb, nzb)
        call Grid_getDeltas(blockList(lb),delta)
        call Grid_getBlkBoundBox(blockList(lb),bndBox)

        delx = delta(IAXIS)
        dely = delta(JAXIS)
        delz = delta(KAXIS)

        delxx = delx * mpole_subSampleInv
        delyy = dely * mpole_subSampleInv
        delzz = delz * mpole_subSampleInv

        ! Get pointer to solution data.

        call Grid_getBlkPtr(blockList(lb), solnData)

        !get limits of blk.. this is compatible with blk sizes not fixed at compile time
        call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

        kmax = blkLimits(HIGH,KAXIS)
        kmin = blkLimits(LOW,KAXIS)  
        jmax = blkLimits(HIGH,JAXIS)
        jmin = blkLimits(LOW,JAXIS)
        imax = blkLimits(HIGH,IAXIS)
        imin = blkLimits(LOW,IAXIS)



        !               Compute potential.

        do k = kmin, kmax
           z = (bndBox(LOW,KAXIS) + (k-kmin)*delz - Zcm) * K3D
           do j = jmin, jmax
              y = (bndBox(LOW,JAXIS) + (j-jmin)*dely - Ycm) * K2D
              do i = imin, imax
                 x = bndBox(LOW,IAXIS) + (i-imin)*delx - Xcm

                 potsum = 0.
                 dvol   = 0.

                 do kk = 1, 1+(mpole_subSample-1)*K3D
                    zz = (z + (kk-0.5)*delzz) * K3D
                    do jj = 1, 1+(mpole_subSample-1)*K2D
                       yy = (y + (jj-0.5)*delyy) * K2D
                       do ii = 1, mpole_subSample
                          xx = x + (ii-0.5)*delxx

                          select case (mpole_geometry)

                          case (G_3DCARTESIAN)
                             ddvol = delxx * delyy * delzz
                             call gr_zonePotential (xx, yy, zz, potential)

                          case (G_3DAXISYMMETRIC)
                             ddvol = delxx * delyy * delzz
                             call gr_zonePotential (xx, yy, zz, potential)

                          case (G_2DCYLINDRICAL)
                             ddvol = twopi * xx * delxx * delyy
                             call gr_zonePotential (xx, 0., yy, potential)
                             potential = potential + point_mass/sqrt(xx**2+yy**2)
                          case (G_1DSPHERICAL)
                             ddvol = fourpi * xx**2 * delxx
                             call gr_zonePotential (xx, 0., 0., potential)
                             potential = potential + point_mass/xx
                          end select

                          potsum = potsum + potential*ddvol
                          dvol   = dvol + ddvol

                       enddo
                    enddo
                 enddo

                 solnData(ipotvar,i,j,k) = -mpfactor*potsum/dvol

              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockList(lb), solnData)

     end do
     
  else
     
     ! 2-D spherical solver
     
     Moment(:,:,1,:,:) = 0.
     
     do lb = 1, blockCount
        call gr_mpoleSphABTerms(blockList(lb))
     end do
     
     ! sum over global j-index

     do m = 0, mpole_mmax
        do l = m, mpole_lmax
           do i = Even, Odd
              call mpi_allreduce (Moment(0,i,1,l,m), Momtmp(0), qmax+1, & 
                   FLASH_REAL, MPI_SUM, gr_meshComm, &
                   error)
              Moment(:,i,1,l,m) = Momtmp(:)
           enddo
        enddo
     enddo
     
     ! sum the innermost contribution over global j-index
     
     call mpi_allreduce (gbnd, gbnd1, 1, & 
          FLASH_REAL, MPI_SUM, gr_meshComm, error)
     gbnd = gbnd1
     
     do l = 0, mpole_lmax, lstep
        
        !      qmax is actually NXMAX+1, The arrays are allocated from 0-qmax
        !      so that when dealing with recurrence relations, the ends don't
        !      have to be treated as special cases for do loop indices.
        !----  but we are really interested in 1:qmax-1 values only
        
        Moment(0,1,1,l,0) = 0.e0       ! this initialization is superfluous
        Moment(qmax-1,2,1,l,0) = 0.e0  ! this is needed and was a bug
        
        do i = 2,qmax                  ! here starting from 2 is right
           Moment(i,1,1,l,0) = &
                &Moment(i,1,1,l,0)+Moment(i,1,2,l,0)*Moment(i-1,1,1,l,0)
        end do
        
        ! ------Here was another bug, I should have started from qmax-2 
        !       and ended at 1
        !       this was may be because of flipping the calculation over
        
        do i = qmax-2,1,-1
           Moment(i,2,1,l,0) = &
                &Moment(i,2,1,l,0)+Moment(i,2,2,l,0)*Moment(i+1,2,1,l,0)
        end do
     end do
     
     ! compute potential for leaf blocks
     
     do lb = 1, blockCount
        call gr_mpoleSphBlkPotential(blockList(lb))
     end do
     
  end if
  
  !===================================================
  
  return
end subroutine gr_mpolePotential

