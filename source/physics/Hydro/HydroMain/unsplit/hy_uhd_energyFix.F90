!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_energyFix
!!
!! NAME
!!
!!  hy_uhd_energyFix
!!
!! SYNOPSIS
!!
!!  call hy_uhd_energyFix( integer (IN) :: blockID,
!!                    integer (IN) :: blkLimits(2,MDIM),
!!                    real(IN)     :: dt,
!!                    real(IN)     :: dtOld,
!!                    real(IN)     :: del(MDIM),
!!                    integer(IN)  :: eosMode)
!!
!! DESCRIPTION
!!  (1) For Hydro and in general:
!!  This routine corrects energy by using the internal energy evolution
!!  to avoid any negativity states of pressure in the Eos routines.
!!
!! (2) For MHD:
!!  This routine corrects energy in two different ways:
!!  The first choice is to fix the energy due to the differences of
!!  the magnetic pressures using the cell-centered magnetic 
!!  fields and the divergence-free cell face-centered magnetic fields.
!!  This correction is optional, but may be useful for low beta
!!  plasma flows. To enable this first correction during the simulation
!!  two runtime parameters "hy_killdivb" and "hy_energyFixSwitch" 
!!  should be both turned on in flash.par file.
!!  The second correction is to use the internal energy evolution
!!  to avoid any negativity states of pressure in the Eos routines.
!!
!!  Note that this routine does more than just correcting energy, and
!!  also computes several quantities such as divergence of magnetic
!!  fields, total pressure, current density, and electric fields, etc.
!!
!! ARGUMENTS
!!
!!  blockID   - a local block ID
!!  blkLimits - an array that holds the lower and upper indices of the section
!!              of block without the guard cells
!!  dt        - time step at n-step
!!  dtOld     - time step at (n-1) step
!!  del       - grid deltas in each direction
!!  eosMode   - a mode used in a call to Eos
!!
!! NOTES
!!
!!  The ENER_VAR and EINT_VAR components of UNK are now always in specific
!!  (i.e., energy per mass) form, on entry as well as on return from this routine.
!!  This is changed from the behavior in FLASH4.2.2 and earlier.
!!
!!  In the MHD case (more specifically when FLASH_UHD_HYDRO is not defined),
!!  the ENER_VAR component of UNK should include magnetic field energy on entry,
!!  but it will not be included on return.
!!***

!!REORDER(4):U


#include "Flash.h"
#include "constants.h"
#include "UHD.h"

Subroutine hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,eosMode)

  use Hydro_data,     ONLY : hy_eswitch, hy_irenorm, hy_geometry,&
                             hy_dtmin, hy_dtminloc, hy_dtminValid, hy_dtminCfl, hy_meshMe, &
                             hy_useAuxEintEqn, hy_smallE,     &
                             hy_cfl, hy_hydroComputeDtOption
#ifdef FLASH_USM_MHD
  use Hydro_data,     ONLY : hy_killdivb, hy_energyFixSwitch
#endif
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_getCellCoords

  implicit none

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt,dtOld
  real, dimension(MDIM), intent(IN) :: del
  integer, intent(IN) :: eosMode
  !! -----------------------------------------------------

  integer :: i,j,k
  real    :: ekin,eint, emag, mhdEnergyCorrection
  real, pointer, dimension(:,:,:,:) :: U
#if !defined FLASH_UHD_HYDRO || !defined FLASH_UHD_3T
    real    :: IntEner,newEint
#endif

  ! Time step calculation
  real :: sndspd2, dt_ltemp, dy, dz
  real, allocatable, dimension(:) :: xCtr, yCtr

  ! cylindrical geometry
  real :: rc, rm, rp

  ! For MHD only
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* For MHD */
  real :: Bxp,Bxm,Byp,Bym,Bzp,Bzm
  real :: cfx2,cfy2,cfz2,bbx2,bby2,bbz2,b2
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: E,Bx,By,Bz
#endif
#endif
#endif /* For MHD */


  !! Get block pointers
  call Grid_getBlkPtr(blockID,U,CENTER)
#ifdef FLASH_USM_MHD /* For MHD */
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_getBlkPtr(blockID,E,SCRATCH)
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif /* For MHD */

  ! initialize geometric terms
  rc = 1.
  rm = 1.
  rp = 1.

  if (hy_geometry /= CARTESIAN) then
     allocate(xCtr(blkLimits(LOW,IAXIS):blkLimits(HIGH, IAXIS)))
     allocate(yCtr(blkLimits(LOW,JAXIS):blkLimits(HIGH, JAXIS)))
     call Grid_getCellCoords&
          (IAXIS,blockID, CENTER,.false.,xCtr, blkLimits(HIGH, IAXIS)-blkLimits(LOW,IAXIS)+1)
     call Grid_getCellCoords&
          (JAXIS,blockId, CENTER,.false.,yCtr, blkLimits(HIGH, JAXIS)-blkLimits(LOW,JAXIS)+1)
  endif


#ifdef FLASH_USM_MHD /* Update cell-centered magnetic fields for USM-MHD */
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           if (hy_energyFixSwitch) then
              !! Store old magnetic pressure from the Godunov update
              U(MAGP_VAR,i,j,k) = &
                   .5*dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),U(MAGX_VAR:MAGZ_VAR,i,j,k))
           endif
#if NDIM > 1
           if (hy_geometry == CYLINDRICAL) then
              rc = abs(xCtr(i))
              rc = 1./rc
              rp = abs(xCtr(i) + 0.5*del(DIR_X))
              rm = abs(xCtr(i) - 0.5*del(DIR_X))
           endif

           U(MAGX_VAR,i,j,k) = &
                .5*(rm*Bx(MAG_FACE_VAR,i,j,k) + rp*Bx(MAG_FACE_VAR,i+1,j,  k  ))*rc

           U(MAGY_VAR,i,j,k) = &
                .5*(   By(MAG_FACE_VAR,i,j,k) +    By(MAG_FACE_VAR,i,  j+1,k  ))
#if NDIM == 3
           U(MAGZ_VAR,i,j,k) = &
                .5*(   Bz(MAG_FACE_VAR,i,j,k) +    Bz(MAG_FACE_VAR,i,  j,  k+1))
#endif
#endif
        enddo
     enddo
  enddo
#endif




  !! Begin do-loops for calculating the followings:
  !! (1) takes care of internal energy,
  !! (2)  - removed -
  !! (3) compute hydro advection dt, and
  !! (4) calculating divergence of velocity fields if requested.
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) .LE. 0.0) then
#endif

              ! In case 3T is used then the energy updates are already done in hy_uhd_unsplitUpdate
              ! and are not needed here.
              ekin = 0.5*dot_product(U(VELX_VAR:VELZ_VAR,i,j,k),U(VELX_VAR:VELZ_VAR,i,j,k))
              emag = 0.0


#ifdef FLASH_USM_MHD /* additional consideration for MHD */
              !! emag is a magnetic pressure calculated using divergence-free
              !! face-centered B fields on the staggered grid and take 
              !! an arithmetic average to get cell-centered B fields 
              !! (see hy_uhd_staggeredDivb.F90)

              !! MAGP_VAR is a magnetic pressure calculated using a standard Godunov
              !! update based on a cell-centered B fields, which are non-divergence-free.
              emag = .5*dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),U(MAGX_VAR:MAGZ_VAR,i,j,k))

              if (hy_killdivb .and. hy_energyFixSwitch) then
                 !! U(MAGP_VAR,i,j,k) is a magnetic pressure from the Godunov cell-centered update
                 mhdEnergyCorrection = emag - U(MAGP_VAR,i,j,k)
                 U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) + mhdEnergyCorrection/U(DENS_VAR,i,j,k)
              endif

              U(MAGP_VAR,i,j,k) = emag
              emag              = emag/U(DENS_VAR,i,j,k)

#endif /* for USM-MHD */

#if !defined FLASH_UHD_HYDRO || !defined FLASH_UHD_3T
              eint = U(ENER_VAR,i,j,k)-ekin-emag

#ifdef EINT_VAR
              if (hy_useAuxEintEqn) then
                 IntEner = U(EINT_VAR,i,j,k)
              end if
#endif

              !! (1) if hy_eswitch -> 0, then we use the total conservative energy to get eint;
              !! (2) if hy_eswitch  > 0, then the eint is obtained by solving a non-conservative
              !!     auxiliary internal energy equation when ener-ekin <= hy_eswitch*(ekin+emag)
              if (.not. hy_useAuxEintEqn .or. eint > hy_eswitch*(ekin+emag)) then
                 !! Specific internal energy
                 newEint = max(hy_smallE,eint)
              else ! Use auxiliary internal energy equation
                 !! Specific internal energy
                 newEint = max(hy_smallE,IntEner)
              endif
              ! Update ener = ekin + eint; note that magp is not included here
              U(ENER_VAR,i,j,k) = newEint + ekin
#ifdef EINT_VAR
              U(EINT_VAR,i,j,k) = newEint
#endif

#endif /* end of #ifndef FLASH_UHD_HYDRO || ifndef FLASH_UHD_3T */



              if (hy_hydroComputeDtOption == 0) then
                 ! ------------------------------------------------------------------------!
                 ! Time step calculation for the next run                                  !
                 ! ------------------------------------------------------------------------!
                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /*compute additional magneto-acoustic speeds for MHD */
                 bbx2 = U(MAGX_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bby2 = U(MAGY_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bbz2 = U(MAGZ_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 b2   = bbx2 + bby2 + bbz2
                 cfx2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbx2))
                 cfy2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bby2))
                 cfz2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbz2))
                 sndspd2 = cfx2
#endif
                 dt_ltemp = hy_cfl*del(DIR_X)/(abs(U(VELX_VAR,i,j,k))+sqrt(sndspd2))
#if NDIM > 1
                 if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then 
                    ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                    dy = del(DIR_Y)
                 else ! Angular coordinates in 2D: Spherical or Polar
                    dy = xCtr(i)*del(DIR_Y)
                 endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sndspd2 = cfy2
#endif
                 dt_ltemp = min(dt_ltemp,hy_cfl*dy/(abs(U(VELY_VAR,i,j,k))+sqrt(sndspd2)))
#if NDIM > 2
                 if (hy_geometry == CARTESIAN) then 
                    ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                    dz = del(DIR_Z)
                 elseif (hy_geometry == CYLINDRICAL) then
                    dz = xCtr(i)*del(DIR_Z)
                 else ! Angular coordinates in 2D: Spherical or Polar
                    dz = xCtr(i)*sin(yCtr(j))*del(DIR_Z) ! z is phi
                 endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sndspd2 = cfz2
#endif
                 dt_ltemp = min(dt_ltemp,hy_cfl*dz/(abs(U(VELZ_VAR,i,j,k))+sqrt(sndspd2)))
#endif /* NDIM > 2 */
#endif /* NDIM > 1 */
                 if( dt_ltemp < hy_dtmin) then ! .or. dt_ltemp < 1.1*dtOld) then
                    hy_dtmin = dt_ltemp
                    hy_dtminloc(1) = i
                    hy_dtminloc(2) = j
                    hy_dtminloc(3) = k
                    hy_dtminloc(4) = blockID
                    hy_dtminloc(5) = hy_meshMe
                    hy_dtminCfl    = hy_cfl
                    hy_dtminValid = .TRUE.
                 end if
                 ! ------------------------------------------------------------------------!
                 ! End of time step calculation for the next run                           !
                 ! ------------------------------------------------------------------------!
              endif


#ifdef DIVV_VAR 
              ! if defined divergence of velocity
              if (hy_geometry == CYLINDRICAL) then
                 rc = abs(xCtr(i))
                 rc = 1./rc
                 rm = abs(xCtr(i-1))
                 rp = abs(xCtr(i+1))
              endif
              U(DIVV_VAR,i,j,k) = &
                   (rp*U(VELX_VAR,i+1,j,k)-rm*U(VELX_VAR,i-1,j,k))/del(DIR_X)*rc
              if (NDIM > 1) then
                 U(DIVV_VAR,i,j,k) = U(DIVV_VAR,i,j,k) &
                      +(U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k))/del(DIR_Y)
                 if (NDIM > 2) then
                    U(DIVV_VAR,i,j,k) = U(DIVV_VAR,i,j,k) &
                         +(U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1))/del(DIR_Z)
                 endif
              endif
              U(DIVV_VAR,i,j,k) = 0.5*U(DIVV_VAR,i,j,k)
#endif /* DIVV_VAR */


              ! Begin additional MHD consideration -------------
#ifndef FLASH_UHD_HYDRO /* For MHD */

           U(DIVB_VAR,i,j,k) =  0.

#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
           if (hy_killdivb) then
              if (NDIM > 1) then
                 Bxp=Bx(MAG_FACE_VAR,i+1,j,    k  )
                 Bxm=Bx(MAG_FACE_VAR,i,  j,    k  )
                 Byp=By(MAG_FACE_VAR,i,  j+1,  k  )
                 Bym=By(MAG_FACE_VAR,i,  j,    k  )
              if (NDIM == 3) then
                 Bzp=Bz(MAG_FACE_VAR,i,  j,    k+1)
                 Bzm=Bz(MAG_FACE_VAR,i,  j,    k  )
              endif
              endif
           else
              U(DIVB_VAR,i,j,k) = 0.
           endif
#endif
#endif
#ifdef FLASH_UGLM_MHD
              if (NDIM > 1) then
                 Bxp=2.0*U(MAGX_VAR,i+1,j,    k  )
                 Bxm=2.0*U(MAGX_VAR,i-1,j,    k  )
                 Byp=2.0*U(MAGY_VAR,i,  j+1,  k  )
                 Bym=2.0*U(MAGY_VAR,i,  j-1,  k  )
              if (NDIM == 3) then
                 Bzp=2.0*U(MAGZ_VAR,i,  j,    k+1)
                 Bzm=2.0*U(MAGZ_VAR,i,  j,    k-1)
              endif
              endif
#endif
           if (hy_geometry == CYLINDRICAL) then
              rc = abs(xCtr(i))
              rc = 1./rc
              rm = abs(xCtr(i) - 0.5*del(DIR_X))
              rp = abs(xCtr(i) + 0.5*del(DIR_X))
           endif
           if (NDIM > 1 .and. hy_killdivb) then
              U(DIVB_VAR,i,j,k) = &
                   (rp*Bxp - rm*Bxm)/del(DIR_X)*rc + (Byp - Bym)/del(DIR_Y)
              if (NDIM == 3) then
                 U(DIVB_VAR,i,j,k) = U(DIVB_VAR,i,j,k) + (Bzp - Bzm)/del(DIR_Z)
              endif
           endif


#ifdef TOTP_VAR
           !! Total plasma pressure and plasma beta
           U(TOTP_VAR,i,j,k) = U(PRES_VAR,i,j,k)+U(MAGP_VAR,i,j,k)
#endif /* TOTP_VAR */
#ifdef BETA_VAR
           U(BETA_VAR,i,j,k) = U(PRES_VAR,i,j,k)/U(MAGP_VAR,i,j,k)
#endif /* BETA_VAR */

#if NFACE_VARS > 0
#if NDIM > 1
#ifdef VECZ_VAR
           !! advance Az to next time step in 2D
           !! dA/dt = -E = VxB - magVisc*J
           U(VECZ_VAR,i,j,k) = U(VECZ_VAR,i,j,k)&
                -0.25*dt*( E(EZ_SCRATCH_GRID_VAR, i,  j,  k) &
                          +E(EZ_SCRATCH_GRID_VAR, i+1,j,  k) &
                          +E(EZ_SCRATCH_GRID_VAR, i,  j+1,k) &
                          +E(EZ_SCRATCH_GRID_VAR, i+1,j+1,k))
#endif /* VECZ_VAR */
#ifdef ELEX_VAR
           U(ELEX_VAR,i,j,k) =( E(EX_SCRATCH_GRID_VAR, i,j,  k  ) &
                               +E(EX_SCRATCH_GRID_VAR, i,j,  k+1) &
                               +E(EX_SCRATCH_GRID_VAR, i,j+1,k  ) &
                               +E(EX_SCRATCH_GRID_VAR, i,j+1,k+1))*0.25
#endif /* ELEX_VAR */
#ifdef ELEY_VAR
           U(ELEY_VAR,i,j,k) =( E(EY_SCRATCH_GRID_VAR, i,  j,k  ) &
                               +E(EY_SCRATCH_GRID_VAR, i+1,j,k  ) &
                               +E(EY_SCRATCH_GRID_VAR, i,  j,k+1) &
                               +E(EY_SCRATCH_GRID_VAR, i+1,j,k+1))*0.25
#endif /* ELEY_VAR */
#ifdef ELEZ_VAR
           U(ELEZ_VAR,i,j,k) =( E(EZ_SCRATCH_GRID_VAR, i,  j,  k) &
                               +E(EZ_SCRATCH_GRID_VAR, i+1,j,  k) &
                               +E(EZ_SCRATCH_GRID_VAR, i,  j+1,k) &
                               +E(EZ_SCRATCH_GRID_VAR, i+1,j+1,k))*0.25
#endif /* ELEZ_VAR */
#endif /* NDIM > 1 */
#endif /* NFACE_VARS > 0 */

#endif /* ifndef FLASH_UHD_HYDRO */
           ! End of additional MHD consideration -------------



#ifdef BDRY_VAR
           endif
#endif

        enddo
     enddo
  enddo


  if (eosMode==MODE_DENS_PRES) then
     U(PRES_VAR,:,:,:) = U(EINT_VAR,:,:,:)*U(DENS_VAR,:,:,:)*(U(GAME_VAR,:,:,:)-1.)
  endif


  !! Release block pointers
  call Grid_releaseBlkPtr(blockID,U,CENTER)
#ifdef FLASH_USM_MHD /* For MHD */
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_releaseBlkPtr(blockID,E,SCRATCH)
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif /* For MHD */


  if (hy_geometry /= CARTESIAN) then
     deallocate(xCtr)
     deallocate(yCtr)
  endif
  return
End Subroutine hy_uhd_energyFix
