!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_redshiftHydro
!!
!!
!! NAME
!!
!!  Cosmology_redshiftHydro
!!
!! SYNOPSIS
!!
!!  Cosmology_redshiftHydro( integer(IN) :: blockCount, 
!!                           integer(IN) :: blockList(blockCount)) 
!!
!! DESCRIPTION
!!
!!  Description:  Applies the redshift operator to hydrodynamical quantities.
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList -   array holding local IDs of blocks on which to advance
!!
!!
!!***

!!REORDER(4): solnData

subroutine Cosmology_redshiftHydro ( blockCount, blockList)
  use Cosmology_data, ONLY :  csm_smlrho, csm_smalle, csm_smallp, &
                              csm_eintSwitch,csm_eintSwitchExist,&
                               csm_oldscaleFactor, csm_scaleFactor, &
             csm_computeRedshiftOnly

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_getBlkIndexLimits
  use Eos_interface, ONLY : Eos_wrapped
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer, INTENT(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList



  real :: gamma, delta_ie, delta_ke
  real :: expfac, oldEner,koldEner, newEner, knewEner
  integer       :: i, j, k,lb
  
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnData(:,:,:,:)

#ifndef PRES_VAR  
! There are no hydro variables to evolve, only particles
  return
#endif

  expfac   = csm_oldscaleFactor / csm_scaleFactor

!------------------------------------------------------------------------

#ifdef PRES_VAR

  ! If we only want the cosmology machinery to tell us what the redshift is, 
  ! exit here and don't do any of this

  if (csm_computeRedshiftOnly) return

  do lb = 1, blockCount
     call Grid_getBlkPtr(blockList(lb),solnData,CENTER)
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              
              gamma = solnData(GAME_VAR,i,j,k)
              koldEner = 0.5 * (solnData(VELX_VAR,i,j,k)**2 + &
                   solnData(VELY_VAR,i,j,k)**2 +&
                   solnData(VELZ_VAR,i,j,k)**2)
              oldEner = solnData(EINT_VAR,i,j,k)
              if ((csm_eintSwitchExist) .and. (csm_eintSwitch >= 0.) .and. &
                   (oldEner >= csm_eintSwitch*koldEner)) &
                   oldEner = solnData(ENER_VAR,i,j,k) - koldEner
              oldEner = max(oldEner, csm_smalle)
              knewEner = koldEner * expfac**4
              newEner = oldEner * expfac**(3.*gamma-1.)
              delta_ie = newEner - oldEner
              delta_ke = knewEner - koldEner
              
              solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) * expfac**2
              solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) * expfac**2
              solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) * expfac**2
              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + &
                                         delta_ie + delta_ke
              solnData(ENER_VAR,i,j,k) = max( solnData(ENER_VAR,i,j,k), csm_smalle )
              solnData(EINT_VAR,i,j,k) = newEner
              solnData(EINT_VAR,i,j,k) = max( solnData(EINT_VAR,i,j,k), csm_smalle )
              solnData(PRES_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) * &
                                         solnData(DENS_VAR,i,j,k) * (gamma-1.)
              solnData(PRES_VAR,i,j,k) = max( solnData(PRES_VAR,i,j,k), csm_smallp )
              
           enddo
        enddo
     enddo
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockList(lb))
     call Grid_releaseBlkPtr(blockList(lb),solnData,CENTER)
     
  enddo

#endif

!-----------------------------------------------------------------------------

  return
end subroutine Cosmology_redshiftHydro

