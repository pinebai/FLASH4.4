!!****if* source/physics/Hydro/HydroMain/unsplit/multiTemp/hy_uhd_multiTempAfter
!!
!! NAME
!!
!!  hy_uhd_multiTempAfter
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_multiTempAfter(integer, intent(in) :: blockCount,
!!                             integer, intent(in) :: blockList(blockCount),
!!                             real,    intent(in) :: dt)
!!
!! DESCRIPTION
!! 
!!
!! ARGUMENTS
!!
!!  blockCount      - Block count
!!  blockList       - Block list
!!  dt              - timestep
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine hy_uhd_multiTempAfter(blockCount, blockList, dt)

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Grid_interface, ONLY: Grid_fillGuardCells
  use Grid_interface, ONLY: Grid_getCellCoords
  use Grid_interface, ONLY: Grid_getDeltas

  use Eos_interface, ONLY: Eos_wrapped
  use hy_uhd_interface, ONLY : hy_uhd_ragelike

  use Hydro_data, ONLY: hy_eosModeAfter
  use Hydro_data, ONLY: hy_3TMode
#ifdef FLASH_USM_MHD
  use Hydro_data,      ONLY : hy_hallVelocity
#ifdef FLASH_UHD_3T
  use hy_uhd_interface,ONLY : hy_uhd_getCurrents
  use hy_memInterface, ONLY : hy_memGetBlkPtr,         &
                              hy_memReleaseBlkPtr
#endif
#endif
  
  
  implicit none

  ! Arguments:
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(blockCount)
  real,    intent(in) :: dt

  ! Local Variables:
  integer :: lb, i,j,k, blockID

  integer :: blkLimitsGC(2,MDIM)
  integer :: blkLimits(2,MDIM)

  real, POINTER, DIMENSION(:,:,:,:) :: U
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr 

  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, allocatable :: zcent(:)

  real :: del(MDIM)

  real :: uele_adv
  real :: uion_adv
  real :: urad_adv
  real :: utot_adv
  real :: pele_adv
  real :: pion_adv
  real :: prad_adv
  real :: dens_new
  real :: utot_new
  real :: uele_new
  real :: uion_new
  real :: urad_new
  real :: uele_ragelike
  real :: uion_ragelike
  real :: urad_ragelike

  real :: divv, divve
  real :: pele
  real :: prad

  real, allocatable :: velx(:,:,:)
  real, allocatable :: vely(:,:,:)
  real, allocatable :: velz(:,:,:)

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
  real, allocatable :: velx_e(:,:,:)
  real, allocatable :: vely_e(:,:,:)
  real, allocatable :: velz_e(:,:,:)
  real, allocatable :: Jp(:,:,:,:)
  real, allocatable :: Jm(:,:,:,:)
#endif              
#endif               
  integer, dimension(MDIM) :: datasize

  integer :: range_switch=UPDATE_ALL

  if(hy_3TMode == HY3T_RAGELIKE) return
  if(hy_3TMode == HY3T_ENTROPY) return

  datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

  ! Need guard cell fill to have access to updated velocities...
  call Grid_fillGuardCells(CENTER,ALLDIR)

  do lb = 1, blockCount
     blockID = blockList(lb)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID, U)
     call Grid_getDeltas(blockID, del)

     allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
     allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
     allocate(zcent(blkLimitsGC(HIGH, KAXIS)))

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
   
     if (hy_hallVelocity) then
        allocate(Jp(3,&
          blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
          blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
          blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
        allocate(Jm(3,&
          blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
          blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
          blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
     
        Jp = 0.0
        Jm = 0.0
     end if

#endif
#endif
     
     call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, blkLimitsGC(HIGH, IAXIS))
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, blkLimitsGC(HIGH, JAXIS))
     call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcent, blkLimitsGC(HIGH, KAXIS))

     allocate( velx( &
          blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
          blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
          blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
     allocate( vely( &
          blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
          blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
          blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
     allocate( velz( &
          blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
          blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
          blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )

     velx = 0.5 * ( U(VOLX_VAR,:,:,:) + U(VELX_VAR,:,:,:) )
     vely = 0.5 * ( U(VOLY_VAR,:,:,:) + U(VELY_VAR,:,:,:) )
     velz = 0.5 * ( U(VOLZ_VAR,:,:,:) + U(VELZ_VAR,:,:,:) )

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
     
     if (hy_hallVelocity) then
          
       allocate( velx_e( &
            blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
            blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
            blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
       allocate( vely_e( &
            blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
            blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
            blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )
       allocate( velz_e( &
            blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1, &
            blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1, &
            blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1) )

      ! get current corrections 
       call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
       call hy_uhd_getCurrents(blockID, range_switch, blkLimits,datasize, del, Jp, Jm, 3,&
                               scrch_Ptr)
       call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
      ! apply corrections to vel and store in vel_e
       velx_e(:,:,:) = velx(:,:,:) - 0.5*(jp(1,:,:,:) + jm(1,:,:,:))   
       vely_e(:,:,:) = vely(:,:,:) - 0.5*(jp(2,:,:,:) + jm(2,:,:,:))  
       velz_e(:,:,:) = velz(:,:,:) - 0.5*(jp(3,:,:,:) + jm(3,:,:,:))
     endif 
#endif
#endif
     
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              dens_new = U(DENS_VAR,i,j,k)
              utot_new = U(EINT_VAR,i,j,k) * dens_new
              uele_adv = U(EELE_VAR,i,j,k) * dens_new
              uion_adv = U(EION_VAR,i,j,k) * dens_new
              urad_adv = U(ERAD_VAR,i,j,k) * dens_new
              utot_adv = uele_adv + uion_adv + urad_adv

              call hy_uhd_getPressure(  &
                   uele_adv / dens_new, &
                   uion_adv / dens_new, &
                   urad_adv / dens_new, &
                   U(:,i,j,k), &
                   pele_adv, &
                   pion_adv, &
                   prad_adv)

              ! pele and prad are used in the work term:
              pele = 0.5 * (pele_adv + U(PELE_VAR,i,j,k))
              prad = 0.5 * (prad_adv + U(PRAD_VAR,i,j,k))

              call compute_divv(divv)
              uele_new = uele_adv - dt * pele * divv
              urad_new = urad_adv - dt * prad * divv
              uion_new = utot_new - uele_new - urad_new

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
              if (hy_hallVelocity) then

                 call compute_divve(divve)
                 uele_new = uele_adv - dt * pele * divve
              endif
#endif              
#endif               
              ! Diagnostic variables...FUNI, FUNE, FUNR are used to
              ! indicated where in the domain internal energies went
              ! negative.
#ifdef FUNI_VAR
              U(FUNI_VAR,i,j,k) = 0.0
              if(uion_new < 0.0) U(FUNI_VAR,i,j,k) = 1.0
#endif
#ifdef FUNE_VAR
              U(FUNE_VAR,i,j,k) = 0.0
              if(uele_new < 0.0) U(FUNE_VAR,i,j,k) = 1.0
#endif
#ifdef FUNR_VAR
              U(FUNR_VAR,i,j,k) = 0.0
              if(urad_new < 0.0) U(FUNR_VAR,i,j,k) = 1.0
#endif

              ! CRASH-like solution has been computed. Now compare it
              ! to the RAGE-like solution. We use the RAGE-like
              ! solution as a fall-back.

              call hy_uhd_ragelike(utot_new-utot_adv, U(:,i,j,k), dens_new, &
                   0.0, 0.0, 0.0, &
                   xcent(i), ycent(j), zcent(k), &
                   uele_ragelike, uion_ragelike, urad_ragelike)

              if(uele_new < 0.0 .or. uion_new < uion_ragelike .or. urad_new < 0.0 ) then
                 uele_new = uele_ragelike
                 uion_new = uion_ragelike
                 urad_new = urad_ragelike
              end if

              U(EELE_VAR,i,j,k) = uele_new / dens_new
              U(EION_VAR,i,j,k) = uion_new / dens_new
              U(ERAD_VAR,i,j,k) = urad_new / dens_new

           enddo
        enddo
     enddo

     call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)

     deallocate(xcent)
     deallocate(ycent)
     deallocate(zcent)
     deallocate(velx)
     deallocate(vely)
     deallocate(velz)
#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
     if (hy_hallVelocity) then
        deallocate(Jm)
        deallocate(Jp)
        deallocate(velx_e)
        deallocate(vely_e)
        deallocate(velz_e)
     end if
#endif
#endif

     call Grid_releaseBlkPtr(blockID, U)
  end do


contains
  
  subroutine compute_divv(divv)
    use Hydro_data, ONLY: hy_geometry
    implicit none

    real, intent(out) :: divv

    divv = 0.0
    if (hy_geometry == CARTESIAN) then
       ! Second order central difference:
       ! divv = (velx(i+1,j,k) - velx(i-1,j,k))/(2*del(IAXIS))

       ! Fourth order, central difference:
       divv = &
            -velx(i+2,j,k) + 8*velx(i+1,j,k) &
            +velx(i-2,j,k) - 8*velx(i-1,j,k)
       divv = divv / (12*del(IAXIS))


#if NDIM >= 2
       ! Second order central difference:
       ! divv = divv + &
       !      (vely(i,j+1,k) - vely(i,j-1,k))/(2*del(JAXIS))

       ! Fourth order, central difference:
       divv = divv + ( &
            -vely(i,j+2,k) + 8*vely(i,j+1,k)     &
            +vely(i,j-2,k) - 8*vely(i,j-1,k) ) / &
            (12*del(JAXIS))
#endif
#if NDIM == 3
       ! Second order central difference:
       ! divv =  divv + &
       !      (velz(i,j,k+1) - velz(i,j,k-1))/(2*del(KAXIS))

       ! Fourth order, central difference:
       divv = divv + ( &
            -velz(i,j,k+2) + 8*velz(i,j,k+1)     &
            +velz(i,j,k-2) - 8*velz(i,j,k-1) ) / &
            (12*del(KAXIS))
#endif

    elseif (hy_geometry == CYLINDRICAL) then
       ! ! Second order central difference:
       ! divv = &
       !      (xcent(i+1)*velx(i+1,j,k) - &
       !      xcent(i-1)*velx(i-1,j,k)) / &
       !      (2*del(IAXIS)*xcent(i))

       ! Fourth order central difference:
       divv = ( &
            -xcent(i+2)*velx(i+2,j,k)    &
            +8*xcent(i+1)*velx(i+1,j,k)  &
            -8*xcent(i-1)*velx(i-1,j,k)  &
            +xcent(i-2)*velx(i-2,j,k) ) / &
            (12*del(IAXIS)*xcent(i))

#if NDIM >= 2
       ! ! Second order central difference:
       ! divv = divv + &
       !      (vely(i,j+1,k) - vely(i,j-1,k))/(2*del(JAXIS))

       ! Fourth order, central difference:
       divv = divv + ( &
            -vely(i,j+2,k) + 8*vely(i,j+1,k) &
            +vely(i,j-2,k) - 8*vely(i,j-1,k) ) / &
            (12*del(JAXIS))

#endif
#if NDIM == 3
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
#endif
    elseif (hy_geometry == SPHERICAL) then
       ! ! Second order central difference:
       ! divv = &
       !      (xcent(i+1)**2*velx(i+1,j,k) - &
       !      xcent(i-1)**2*velx(i-1,j,k)) / &
       !      (2*del(IAXIS)*xcent(i)**2)

       ! Fourth order central difference:
       divv = ( &
            -(xcent(i+2)**2)*velx(i+2,j,k)    &
            +8*(xcent(i+1)**2)*velx(i+1,j,k)  &
            -8*(xcent(i-1)**2)*velx(i-1,j,k)  &
            +(xcent(i-2)**2)*velx(i-2,j,k) ) / &
            (12*del(IAXIS)*xcent(i)**2)

#if NDIM >= 2
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
#endif
       
    else
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
    end if

  end subroutine compute_divv

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
  
  subroutine compute_divve(divve)
    use Hydro_data, ONLY: hy_geometry
    implicit none

    real, intent(out) :: divve

    divve = 0.0
    if (hy_geometry == CARTESIAN) then
       ! Second order central difference:
       ! divv = (velx(i+1,j,k) - velx(i-1,j,k))/(2*del(IAXIS))

       ! Fourth order, central difference:
       divve = &
            -velx_e(i+2,j,k) + 8*velx_e(i+1,j,k) &
            +velx_e(i-2,j,k) - 8*velx_e(i-1,j,k)
       divve = divve / (12*del(IAXIS))


#if NDIM >= 2
       ! Second order central difference:
       ! divv = divv + &
       !      (vely(i,j+1,k) - vely(i,j-1,k))/(2*del(JAXIS))

       ! Fourth order, central difference:
       divve = divve + ( &
            -vely_e(i,j+2,k) + 8*vely_e(i,j+1,k)     &
            +vely_e(i,j-2,k) - 8*vely_e(i,j-1,k) ) / &
            (12*del(JAXIS))
#endif
#if NDIM == 3
       ! Second order central difference:
       ! divv =  divv + &
       !      (velz(i,j,k+1) - velz(i,j,k-1))/(2*del(KAXIS))

       ! Fourth order, central difference:
       divve = divve + ( &
            -velz_e(i,j,k+2) + 8*velz_e(i,j,k+1)     &
            +velz_e(i,j,k-2) - 8*velz_e(i,j,k-1) ) / &
            (12*del(KAXIS))
#endif

    elseif (hy_geometry == CYLINDRICAL) then
       ! ! Second order central difference:
       ! divv = &
       !      (xcent(i+1)*velx(i+1,j,k) - &
       !      xcent(i-1)*velx(i-1,j,k)) / &
       !      (2*del(IAXIS)*xcent(i))

       ! Fourth order central difference:
       divve = ( &
            -xcent(i+2)*velx_e(i+2,j,k)    &
            +8*xcent(i+1)*velx_e(i+1,j,k)  &
            -8*xcent(i-1)*velx_e(i-1,j,k)  &
            +xcent(i-2)*velx_e(i-2,j,k) ) / &
            (12*del(IAXIS)*xcent(i))

#if NDIM >= 2
       ! ! Second order central difference:
       ! divv = divv + &
       !      (vely(i,j+1,k) - vely(i,j-1,k))/(2*del(JAXIS))

       ! Fourth order, central difference:
       divve = divve + ( &
            -vely_e(i,j+2,k) + 8*vely_e(i,j+1,k) &
            +vely_e(i,j-2,k) - 8*vely_e(i,j-1,k) ) / &
            (12*del(JAXIS))

#endif
#if NDIM == 3
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
#endif
    elseif (hy_geometry == SPHERICAL) then
       ! ! Second order central difference:
       ! divv = &
       !      (xcent(i+1)**2*velx(i+1,j,k) - &
       !      xcent(i-1)**2*velx(i-1,j,k)) / &
       !      (2*del(IAXIS)*xcent(i)**2)

       ! Fourth order central difference:
       divve = ( &
            -(xcent(i+2)**2)*velx_e(i+2,j,k)    &
            +8*(xcent(i+1)**2)*velx_e(i+1,j,k)  &
            -8*(xcent(i-1)**2)*velx_e(i-1,j,k)  &
            +(xcent(i-2)**2)*velx_e(i-2,j,k) ) / &
            (12*del(IAXIS)*xcent(i)**2)

#if NDIM >= 2
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
#endif
       
    else
       call Driver_abortFlash("[hy_uhd_getRiemannState] Unsupported geometry for computing DIVV")
    end if

  end subroutine compute_divve
#endif
#endif
  
end subroutine hy_uhd_multiTempAfter
