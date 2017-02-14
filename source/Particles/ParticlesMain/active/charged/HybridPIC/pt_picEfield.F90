!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picEfield
!!
!! NAME
!!
!!  pt_picEfield
!!
!! SYNOPSIS
!!
!!  call pt_picEfield(integer(in) :: ibx,
!!                    integer(in) :: iby,
!!                    integer(in) :: ibz)
!!
!! DESCRIPTION
!!
!!  Compute the electric field from the magnetic field in u(ibxyz,:,:,:)
!!
!! ARGUMENTS
!!
!!   ibx : variable index for Bx
!!
!!   iby : variable index for By
!!
!!   ibz : variable index for Bz
!!
!!
!!
!!***

subroutine pt_picEfield(ibx, iby, ibz)

  use pt_picData, ONLY: pt_picGam, pt_picTe, pt_picResistivityHyper, &
       pt_picCdensMin, pt_picPdensity_1, pt_picPdensity_2, pt_picPcharge_1, &
       pt_picPcharge_2  
  use pt_picInterface, ONLY : pt_picCurl
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,&
       Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits, Grid_releaseBlkPtr, &
       Grid_fillGuardCells

  use Particles_data, ONLY: pt_meshMe

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  implicit none

  integer, intent(in) :: ibx, iby, ibz
  real                :: pec
  integer :: blockList(MAXBLOCKS), blockCount, localSize(MDIM)
  real :: blockSize(MDIM), blockCenter(MDIM)
  real :: je(MDIM), jexb(MDIM), h(MDIM), dv, b(3), cden_inv
  integer :: blkLimits(LOW:HIGH,MDIM), blkLimitsGC(LOW:HIGH,MDIM)
  real, dimension(:,:,:,:), pointer :: u
  integer :: i, j, k, block_no

  real, parameter :: ec  = 1.6021773e-19 ! Elementary charge [C]
  ! permeability of vacuum [Vs/Am] (magnetic constant) 
  real, parameter :: mu0 = 4.0e-7*3.141592653589793238462643383279503
  real, parameter :: kb = 1.38066e-23 ! Boltzmann [J/K = kg(m/s)^2/K]

  real :: cmin             ! minimum charge density
  real :: gres             ! resistivity to use
  cmin = (pt_picPdensity_1*pt_picPcharge_1 &
       +pt_picPdensity_2*pt_picPcharge_2)*ec*pt_picCdensMin

  pec = kb/ec*(ec*pt_picPdensity_1)**(1-pt_picGam)*pt_picTe

  ! compute J = rot(B)/mu0. Assumes that rot(external field) = 0
  call pt_picCurl(ibx, iby, ibz, GRJX_VAR, GRJY_VAR, GRJZ_VAR, 1.0/mu0, .true.)

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  do block_no = 1, blockCount
     call Grid_getBlkPtr(blockList(block_no), u)
     call Grid_getBlkPhysicalSize(blockList(block_no), blockSize)
     call Grid_getBlkIndexLimits(blockList(block_no), &
          blkLimits, blkLimitsGC, CENTER)
     localSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1   ! NXB, NYB and NZB

     h = blockSize/localSize  ! cell size
     dv = product(h(1:NDIM))  ! cell volume [m^3]

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              ! e- current, Je=J-Ji
              je(1) = u(GRJX_VAR,i,j,k)-u(GJIX_VAR,i,j,k)
              je(2) = u(GRJY_VAR,i,j,k)-u(GJIY_VAR,i,j,k)
              je(3) = u(GRJZ_VAR,i,j,k)-u(GJIZ_VAR,i,j,k)
              ! compute Je x B
              b(1) = u(ibx,i,j,k) + u(GBX1_VAR,i,j,k)
              b(2) = u(iby,i,j,k) + u(GBY1_VAR,i,j,k)
              b(3) = u(ibz,i,j,k) + u(GBZ1_VAR,i,j,k)
              jexb(1) = je(2)*b(3) - je(3)*b(2)
              jexb(2) = je(3)*b(1) - je(1)*b(3)
              jexb(3) = je(1)*b(2) - je(2)*b(1)
              ! Add -Grad(pe) = -pec*Grad(cden^gam)
              if (pt_picGam > 0.0) then
                 jexb(1) = jexb(1) - pec*0.5*(u(CDEN_VAR,i+1,j,k)**pt_picGam &
                      -u(CDEN_VAR,i-1,j,k)**pt_picGam)/h(1)
                 if (NDIM > 1) then
                    jexb(2) = jexb(2) - pec*0.5*(u(CDEN_VAR,i,j+1,k)**pt_picGam &
                         -u(CDEN_VAR,i,j-1,k)**pt_picGam)/h(2)
                 end if
                 if (NDIM == 3) then
                    jexb(3) = jexb(3) - pec*0.5*(u(CDEN_VAR,i,j,k+1)**pt_picGam &
                         -u(CDEN_VAR,i,j,k-1)**pt_picGam)/h(3)
                 end if
              end if
              gres = u(GRES_VAR,i,j,k)
              if (u(CDEN_VAR,i,j,k) < cmin) then ! vacuum
                 cden_inv = 1.0/cmin
              else
                 cden_inv = 1.0/u(CDEN_VAR,i,j,k)
              end if
              u(GREX_VAR,i,j,k) = jexb(1)*cden_inv
              u(GREY_VAR,i,j,k) = jexb(2)*cden_inv
              u(GREZ_VAR,i,j,k) = jexb(3)*cden_inv
              ! E-field w/o resistive term for Lorentz force
              u(GEPX_VAR,i,j,k) = u(GREX_VAR,i,j,k)
              u(GEPY_VAR,i,j,k) = u(GREY_VAR,i,j,k)
              u(GEPZ_VAR,i,j,k) = u(GREZ_VAR,i,j,k)
              ! + neta*J for resistivity
              u(GREX_VAR,i,j,k) = u(GREX_VAR,i,j,k) &
                   + gres*u(GRJX_VAR,i,j,k)
              u(GREY_VAR,i,j,k) = u(GREY_VAR,i,j,k) &
                   + gres*u(GRJY_VAR,i,j,k)
              u(GREZ_VAR,i,j,k) = u(GREZ_VAR,i,j,k) &
                   + gres*u(GRJZ_VAR,i,j,k)
              ! - neta2*Lap(J) for hyperresistivity
              if (NDIM == 1) then
                 u(GREX_VAR,i,j,k) = u(GREX_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJX_VAR,i+1,j,k) -2*u(GRJX_VAR,i,j,k) &
                      + u(GRJX_VAR,i-1,j,k))/(h(1)*h(1))
                 u(GREY_VAR,i,j,k) = u(GREY_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJY_VAR,i+1,j,k) -2*u(GRJY_VAR,i,j,k) &
                      + u(GRJY_VAR,i-1,j,k))/(h(1)*h(1))
                 u(GREZ_VAR,i,j,k) = u(GREZ_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJZ_VAR,i+1,j,k) -2*u(GRJZ_VAR,i,j,k) &
                      + u(GRJZ_VAR,i-1,j,k))/(h(1)*h(1))
              end if
              if (NDIM > 1) then
                 u(GREX_VAR,i,j,k) = u(GREX_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJX_VAR,i,j+1,k) -2*u(GRJX_VAR,i,j,k) &
                      + u(GRJX_VAR,i,j-1,k))/(h(2)*h(2))
                 u(GREY_VAR,i,j,k) = u(GREY_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJY_VAR,i,j+1,k) -2*u(GRJY_VAR,i,j,k) &
                      + u(GRJY_VAR,i,j-1,k))/(h(2)*h(2))
                 u(GREZ_VAR,i,j,k) = u(GREZ_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJZ_VAR,i,j+1,k) -2*u(GRJZ_VAR,i,j,k) &
                      + u(GRJZ_VAR,i,j-1,k))/(h(2)*h(2))
              end if
              if (NDIM == 3) then
                 u(GREX_VAR,i,j,k) = u(GREX_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJX_VAR,i,j,k+1) -2*u(GRJX_VAR,i,j,k) &
                      + u(GRJX_VAR,i,j,k-1))/(h(3)*h(3))
                 u(GREY_VAR,i,j,k) = u(GREY_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJY_VAR,i,j,k+1) -2*u(GRJY_VAR,i,j,k) &
                      + u(GRJY_VAR,i,j,k-1))/(h(3)*h(3))
                 u(GREZ_VAR,i,j,k) = u(GREZ_VAR,i,j,k) &
                      - pt_picResistivityHyper* &
                      (u(GRJZ_VAR,i,j,k+1) -2*u(GRJZ_VAR,i,j,k) &
                      + u(GRJZ_VAR,i,j,k-1))/(h(3)*h(3))
              end if
              ! + C not implemented
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(block_no), u)
  end do

  call Grid_fillGuardCells( CENTER, ALLDIR)

end subroutine pt_picEfield
