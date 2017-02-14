!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commCheckRayExchange
!!
!!  NAME     
!!   ed_commCheckRayExchange
!!
!!  SYNOPSIS
!!   ed_commCheckRayExchange()
!!
!!  DESCRIPTION 
!!
!!  ARGUMENTS
!!
!!  NOTES
!!
!!  SIDE EFFECTS
!!
!!***

subroutine ed_commCheckRayExchange()
  use EnergyDeposition_data, ONLY : ed_rays, ed_rayCount
  use ed_commData, ONLY : ed_commLog, ed_commLogUnit
  use UTPipeline, ONLY : UTPipeline_iterateItems, UTPipeline_isEmpty
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  interface
     subroutine ed_commLogRay(ray, rayDescription)
       implicit none
       real, dimension(:), intent(IN) :: ray
       character(len=*), intent(IN) :: rayDescription
     end subroutine ed_commLogRay
  end interface
  integer :: i
  character(len=200) :: miscString
  character(len=*), parameter :: badRayExchangeMsg = &
       'There should not be any active rays! See local log files.'
  logical :: isPipelineEmpty, isBadRayExchange

  call UTPipeline_isEmpty(isPipelineEmpty)
  isBadRayExchange = (ed_rayCount > 0 .or. .not.isPipelineEmpty)

  if (isBadRayExchange) then
     if (.not.ed_commLog) call Logfile_open(ed_commLogUnit,.true.)
     write(ed_commLogUnit,'(a)') 'There are still active rays. This is bad.'
     call UTPipeline_iterateItems(ed_commLogRay)
     do i = 1, ed_rayCount
        write (miscString,'(a,i10)') 'ed_rays: ray ', i
        call ed_commLogRay(ed_rays(:,i), trim(miscString))
     end do
     call Logfile_close(.true.)
     call Driver_abortFlash(badRayExchangeMsg)
  end if
end subroutine ed_commCheckRayExchange

#include "EnergyDeposition.h"
subroutine ed_commLogRay(ray, rayDescription)
  use ed_commData, ONLY : ed_commLogUnit
  implicit none
  real, dimension(:), intent(IN) :: ray
  character(len=*), intent(IN) :: rayDescription
  character(len=*), parameter :: logStr = "(a, 3(/,a,i10), 3(/,a,es25.18))"
  write(ed_commLogUnit,logStr) &
       rayDescription, &
       "  Proc ", int(ray(RAY_PROC)), &
       "  Block", int(ray(RAY_BLCK)), &
       "  Tag  ", int(ray(RAY_TAGS)), &
       "  Posx ", ray(RAY_POSX), &
       "  Posy ", ray(RAY_POSY), &
       "  Posz ", ray(RAY_POSZ)
end subroutine ed_commLogRay
