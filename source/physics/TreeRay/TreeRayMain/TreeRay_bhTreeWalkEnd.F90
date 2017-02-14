!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhTreeWalkEnd
!!
!! NAME
!!
!!  TreeRay_bhTreeWalkEnd
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeRay_bhTreeWalkEnd(iterate)
  use TreeRay_data, ONLY : tr_bhMaxRelEradErr, tr_bhLocRelErr, tr_meshMe, &
    tr_comm, tr_bhRelErr, tr_bhErrControl, tr_bhEradTot, tr_bhLocEradTot, &
    tr_bhOldEradTot, tr_bhMionTot, tr_bhLocMionTot, tr_bhOldMionTot, &
    tr_bhUseTreeRay

  use Logfile_interface, ONLY : Logfile_stamp    
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  logical, intent(OUT) :: iterate
  integer :: ierr
  character(len=MAX_STRING_LENGTH) :: strBuff  
  real :: RelEradErr, RelMionErr

  if (.not. tr_bhUseTreeRay) then
    iterate = .false.
    return
  endif

  select case (tr_bhErrControl)

    case ("erad_cell")
      !call MPI_AllReduce(tr_bhLocRelErr,tr_bhMaxRelEradErr,1, &
      !& FLASH_REAL,MPI_MAX,tr_comm, ierr)
      call MPI_AllReduce(tr_bhLocRelErr,tr_bhMaxRelEradErr,1, &
      & FLASH_REAL,MPI_MAX, tr_comm, ierr)
  
      if (tr_meshMe == MASTER_PE) then
        write (strBuff, '("Relative erad error: ", e12.5)') tr_bhMaxRelEradErr
        call Logfile_stamp( strBuff, "[TreeRay]")
      endif

      if (tr_bhMaxRelEradErr > tr_bhRelErr) then
        iterate = .true.
      else
        iterate = .false.
      endif

    case ("erad_tot")
      tr_bhOldEradTot = tr_bhEradTot
      call MPI_AllReduce(tr_bhLocEradTot,tr_bhEradTot,1, &
      & FLASH_REAL,FLASH_SUM,tr_comm, ierr)
      RelEradErr = 2.0*abs(tr_bhEradTot - tr_bhOldEradTot) &
      &          / (tr_bhEradTot + tr_bhOldEradTot + 1d-99)

      if (tr_meshMe == MASTER_PE) then
        write (strBuff, '("Relative erad error: ", e12.5)') RelEradErr
        call Logfile_stamp( strBuff, "[TreeRay]")
      endif

      if (RelEradErr > tr_bhRelErr) then
        iterate = .true.
      else
        iterate = .false.
      endif

    case ("mion_tot")
      tr_bhOldMionTot = tr_bhMionTot
      call MPI_AllReduce(tr_bhLocMionTot,tr_bhMionTot,1, &
      & FLASH_REAL,FLASH_SUM,tr_comm, ierr)
      RelMionErr = 2.0*abs(tr_bhMionTot - tr_bhOldMionTot) &
      &          / (tr_bhMionTot + tr_bhOldMionTot + 1d-99)

      if (tr_meshMe == MASTER_PE) then
        write (strBuff, '("Relative erad error: ", e12.5)') RelMionErr
        call Logfile_stamp( strBuff, "[TreeRay]")
      endif

      if (RelMionErr > tr_bhRelErr) then
        iterate = .true.
      else
        iterate = .false.
      endif


    case default
      call Driver_abortFlash("Unknow error control methon in TreeRay_bhTreeWalkEnd.")
  end select


  return
end subroutine TreeRay_bhTreeWalkEnd

