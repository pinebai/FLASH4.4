!!****if* source/Grid/GridMain/Grid_renormAbundance
!!
!! NAME
!!
!!  Grid_renormAbundance
!!
!!
!! SYNOPSIS
!!
!!  Grid_renormAbundance(integer(IN) :: blockID,
!!                     integer(IN) :: blkLimits(2,MDIM),
!!                     real,pointer :: solnData(:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Renormalize the abundances in a given block so they sum to 1.
!!  This should be used before calling the EOS.  Each abundance is
!!  restricted to fall between smallx and 1.
!!
!!  Also check the abundance conservation and report if it is worse
!!  than abundErr.
!!
!!  This routine is called automatically by Hydro and MHD in FLASH 
!!  if the irenorm runtime parameter is set to 1.
!!
!!  Only abundances/fluids which contribute to the EOS are included.
!!
!!
!! ARGUMENTS
!!
!!  blockID -   the block number to renormalize
!!  blkLimits - the index limits for internal zones of the block to renormalize
!!  solnData -  Pointer to the block to be renormalized
!!
!!
!! PARAMETERS
!!
!!  smallx -    the cutoff value for the composition mass fraction
!!
!!
!! SEE ALSO 
!!
!!  Grid_limitAbundance
!!
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


#ifdef DEBUG
#define DEBUG_MS
#endif

subroutine Grid_renormAbundance(blockId,blkLimits,solnData)


  use Grid_data, ONLY : gr_smallx
  use Grid_interface, ONLY : Grid_getSingleCellCoords
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(in) :: blockId
  integer, intent(in), dimension(2,MDIM)::blkLimits
  real,pointer :: solnData(:,:,:,:)

  integer :: i, j, k, n
  
  real :: sum, suminv, error

  real, parameter :: abundErr = 1.e-4
  integer,dimension(MDIM) :: point
  real,dimension(MDIM)::pntCoord
  character(len=120)::log_message


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
               
           sum = 0.e0

           ! loop over all of the abundances in 
           !the current zone and retrict them to
           ! fall between smallx and 1.  Then 
           !compute the sum of the abundances

           do n = SPECIES_BEGIN,SPECIES_END

              solnData(n,i,j,k) = & 
                   max(gr_smallx, &
                   min(1.e0,solnData(n,i,j,k)))
              sum = sum + solnData(n,i,j,k)
              
           enddo

! if the conservation is really messed up, give an error
           error = abs(sum - 1.e0)
!           if (error .GT. abundErr) then
           if (error >= 0.10) then
              point(IAXIS)=i; point(JAXIS)=j; point(KAXIS)=k
              print*,'get grid single cell coords',point
              call Grid_getSingleCellCoords&
                   (point,blockId,CENTER, EXTERIOR, pntCoord)
              print *, 'Error: non-conservation in block ', blockId
              print *, 'Abundance non-conservation by ', error
              
              print *, 'x = ', pntCoord(IAXIS)
              print *, 'y = ', pntCoord(JAXIS)
              print *, 'z = ', pntCoord(KAXIS)
#ifdef DENS_VAR
              print *, 'density = ', solnData(DENS_VAR,i,j,k)
#endif
              print *, 'xnuc = ', solnData(SPECIES_BEGIN:&
                   SPECIES_END,i,j,k)

!!$              write(log_message, &
!!$                   '(a, g15.8, a, g15.8, a, g15.8, a, g15.8)')  &
!!$                   '!! non-cons. by ', error, ' at ', pntCoord 
!!$              call Logfile_stamp(log_message)

! bail if the error is exceptionally large -- something is seriously wrong
              if (error > .10) then
                 call Driver_abortFlash('Error too high in abundances')
              endif
              
           endif
           
           ! compute the inverse of the sum and multiply all of the abundances by
           ! this value to get the abundances summing to 1 once again
           suminv = 1.e0 / sum
           
           do n = SPECIES_BEGIN, SPECIES_END
              solnData(n,i,j,k) =  & 
                   max(gr_smallx, min(1.e0,suminv*&
                   solnData(n,i,j,k)))
           enddo
               
        enddo
     enddo
  enddo
  
!==============================================================================
  return
end subroutine Grid_renormAbundance

