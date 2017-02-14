!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Ionize
!!
!! NAME
!!  Ionize
!!
!! SYNOPSIS
!!
!!  call Ionize(integer(IN) :: blockCount,
!!              integer(IN) :: blockList(blockCount),
!!              real(IN)    :: dt,
!!              real(IN)    :: time)
!!
!!
!!
!! DESCRIPTION
!! Apply the ionization operator
!! on the list of blocks provided as input
!!
!! ARGUMENTS
!! blockCount : The number of blocks in the list
!! blockList(:) : The list of blocks on which to apply the stirring operator
!! dt : the current timestep
!! time : the current time
!!
!! NOTES
!!          In the present version there is no evaluation of the
!!          energetic contribution of the ionization and recombination
!!          to the energy equation (heating and cooling the plasma). We
!!          are going to implement this computation soon.
!!
!!          The source term we considered is adequate to solve the
!!          problem for optically thin plasma in the ``coronal''
!!          approximation. This means that we are considering
!!          collisional ionization, auto-ionization, radiative
!!          recombination and dielectronic recombination. The
!!          possibility to change the source term is foreseen to treat
!!          other problems (e.g. the nebular model). We used the
!!          ionization and recombination coefficients computed by
!!          Summers.
!!
!!***

!======================================================================
subroutine Ionize(blockCount,blockList,dt,time)
  use Ionize_data, ONLY : ion_xfrac,ion_dneimax,ion_dneimin,&
       ion_tneimin,ion_tneimax,ion_idx,ion_smallx,ion_ELEC, &
       ion_emass, ion_useIonize
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  integer :: blk, i, j, k, n, nphases, blockID
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC
!
  integer :: nblocks, ierr
  real,pointer, dimension(:,:,:,: ) :: solnData
!
  real, dimension(NSPECIES) :: xin, xout
  !
  real :: sdot
  !
  real :: tmp, rho, ei, ek, den  ! den is also defined in Ionize_data... DEV suspicious
!DEV: CD - I think den is intended as a local variable, so have reinserted den variable.
  logical :: nei_zone

!
!==============================================================================
!
  if (.NOT. ion_useIonize) return !  RETURN immediately

  do blk = 1, blockCount
     blockID = blockList(blk)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     ! get the ratio xfrac = [mass fraction]/[population fraction]
     !
     ! get the current timestep

     nei_zone = .FALSE.
!
! sweep over all the zones
!
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!
              tmp  = solnData(TEMP_VAR,i,j,k)
              rho  = solnData(DENS_VAR,i,j,k)
!
              sdot = 0.0e0
!
! compute the electron number density
!
              den   = rho*ion_xfrac(ion_ELEC)/ion_emass
!
! if the temperature is between the limits and
! the density is between the limits
!
              if ( (tmp >= ion_tneimin .AND. tmp <= ion_tneimax) .AND.         &
                   (den >= ion_dneimin .AND. den <= ion_dneimax) ) then
!
                 nei_zone = .TRUE.
!
! load the population fractions DEV :the start and end points may need fix
                 xin = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
!
! evolve them                  
                 call neimn(dt, tmp, den, xin, xout, sdot, ion_idx, ion_xfrac)
!
! update the global population fraction arrays
                 do n = 1, NSPECIES
                    xout(n) = max(xout(n), ion_smallx)
                 enddo
!
                  solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k) = xout
!
! change in internal energy due to energy released
                  if (sdot.ne.0.0E0) then 
                     ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +             &
     &                           solnData(VELY_VAR,i,j,k)**2 +             &
     &                           solnData(VELZ_VAR,i,j,k)**2)
!
                     ei = solnData(ENER_VAR,i,j,k) - ek
                     ei = ei + dt*sdot
!
! update the global thermodynamic quantities due to the nei evolution
                     solnData(EINT_VAR,i,j,k) = ei
                     solnData(ENER_VAR,i,j,k) = ei + ek
!
#ifdef ENUC_VAR
                     solnData(ENUC_VAR,i,j,k) = sdot
#endif
                  endif
               endif
!
            enddo
         enddo
      enddo

      call Grid_releaseBlkPtr(blockID,solnData)
!
! if we called the NEIMN on any zones in this block, then crank the
! eos out on the entire block
      if (nei_zone) then 
         call Eos_wrapped(MODE_DENS_EI,blkLimitsGC,blockID)
      endif

   end do
   return
 end subroutine Ionize


