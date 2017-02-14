!!****if* source/physics/Hydro/HydroMain/Hydro_recalibrateEints
!! NAME
!!
!!  Hydro_recalibrateEints
!! 
!! SYNOPSIS
!!
!!  call Hydro_recalibrateEints(
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID)
!!
!!  call Hydro_recalibrateEintsForCell(
!!                     integer(in)    :: eint,
!!                     integer(INOUT) :: eion,
!!                     integer(INOUT) :: eele,
!!            OPTIONAL,integer(INOUT) :: erad,
!!            OPTIONAL,integer(INOUT) :: e1,e2,e3)
!!  call Hydro_recalibrateEints(
!!                     integer(in)    :: eint,
!!                     integer(INOUT) :: eion,
!!                     integer(INOUT) :: eele,
!!            OPTIONAL,integer(INOUT) :: erad,
!!            OPTIONAL,integer(INOUT) :: e1,e2,e3)
!!
!!
!! DESCRIPTION
!!
!! This function recalibrates multiTemp component internal energies and energies
!! so that their sums agree with the overall values for eint.
!!
!! In applications not configured for multiTemp Hydro, the call will act like
!! stubs and leave all arguments unchanged.
!!
!! The assumption is that
!!
!!    eint = eion + eele + erad    (in general)
!!    eint = eion + eele           (in the one-cell form if erad is not present;
!!                                  in the per-range/block form if runtime parameter
!!                                  hy_3Ttry_D==3.0).
!!
!! In the per-range/block form of the routine, those variables represent
!! the permanent solution variables (in UNK) EINT_VAR, EION_VAR, EELE_VAR,
!! and ERAD_VAR.
!!
!! If the equation is not true, or is only approximately try, the call will make
!! it true if possible by multiplying the provided two or three component energies
!! with a common scaling factor.
!!
!! Hydro_recalibrateEints is available as a generic function.
!! There is a specific implementation that acts on all cells within a range for
!! a given block.
!! Another specific implementation, available with the specific name
!! Hydro_recalibrateEintsForCell, can be used to adjust energies fo the
!! components of one cell.
!!
!!  ARGUMENTS 
!!
!!   
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   blockID: current block number
!!
!!   eint     : combined internal energy (input)
!!   eion,eele: "ion" and "electron" component energies; to be adjusted.
!!   erad     : radiation component energy; to be adjusted, if present.
!!   e1,e2,e3 : alternative energy variables to be co-adjusted with eion,eele,
!!              erad, respectively,  if present (may not be useful)
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then recalibration is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then recalibration is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Hydro_interface. All functions calling this routine should include
!!      a statement like
!!      use Hydro_interface, ONLY : Hydro_recalibrateEints
!!
!!
!!  SEE ALSO
!!
!!     Hydro_interface
!!     Eos_wrapped
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData
#include "Hydro_components.h"

subroutine Hydro_recalibrateEints(range,blockID)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr

#include "Flash.h"

#if HYDRO_NUM_EINT_COMPONENTS > 1
#ifdef FLASH_UHD_3T
  use hy_uhd_MultiTempData, ONLY : hy_3Ttry_D
#else
  use Hydro_data, ONLY : hy_3Ttry_D
#endif
#define SKIP_ERAD (hy_3Ttry_D == 3.0)
#else
#define SKIP_ERAD (.TRUE.)
#endif

  implicit none

#include "constants.h"

  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID

#if HYDRO_NUM_COMPONENTS > 1
  real, pointer:: solnData(:,:,:,:)

  real :: sumEi, scaleEi, prevEi
  integer :: dataStruct
  integer :: i,j,k

!! ---------------------------------------------------------------------------------


  dataStruct=CENTER
  call Grid_getBlkPtr(blockID,solnData,dataStruct)

  do k = range(LOW,KAXIS), range(HIGH,KAXIS)
     do j = range(LOW,JAXIS), range(HIGH,JAXIS)
        do i = range(LOW,IAXIS), range(HIGH,IAXIS)
           

           sumEi = solnData(EION_VAR,i,j,k)
           if (HYDRO_NUM_EINT_COMPONENTS .GE. 2) sumEi = sumEi + solnData(EELE_VAR,i,j,k)
           if (HYDRO_NUM_EINT_COMPONENTS .GE. 3 .AND. .NOT. SKIP_ERAD) sumEi = sumEi + solnData(ERAD_VAR,i,j,k)
           scaleEi = solnData(EINT_VAR,i,j,k) / sumEi
           

           prevEi = solnData(EION_VAR,i,j,k)
           solnData(EION_VAR,i,j,k) = solnData(EION_VAR,i,j,k) * scaleEi
#ifdef E1_VAR
           solnData(E1_VAR,i,j,k) = solnData(E1_VAR,i,j,k) - prevEi + solnData(EION_VAR,i,j,k)
#endif
           
           if (HYDRO_NUM_EINT_COMPONENTS .GE. 2) then
              prevEi = solnData(EELE_VAR,i,j,k)
              solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) * scaleEi
#ifdef E2_VAR
              if (HYDRO_NUM_E_COMPONENTS .GE. 2) then
                 solnData(E2_VAR,i,j,k) = solnData(E2_VAR,i,j,k) - prevEi + solnData(EELE_VAR,i,j,k)
              end if
#endif
           end if

           if (HYDRO_NUM_EINT_COMPONENTS .GE. 3 .AND. .NOT. SKIP_ERAD) then
              prevEi = solnData(ERAD_VAR,i,j,k)
              solnData(ERAD_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k) * scaleEi
#ifdef E3_VAR
              if (HYDRO_NUM_E_COMPONENTS .GE. 3) then
                 solnData(E3_VAR,i,j,k) = solnData(E3_VAR,i,j,k) - prevEi + solnData(ERAD_VAR,i,j,k)
              end if
#endif
           end if
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,dataStruct)

#endif

  return
end subroutine Hydro_recalibrateEints

subroutine Hydro_recalibrateEintsForCell(eint,eion,eele,erad,e1,e2,e3)
  implicit none
  real,intent(in)    :: eint
  real,intent(INOUT) :: eion,eele
  real,intent(INOUT),OPTIONAL :: erad
  real,intent(INOUT),OPTIONAL :: e1,e2,e3

#if HYDRO_NUM_EINT_COMPONENTS > 1
  real :: sumEi, scaleEi, prevEi

  sumEi = eion
  if (HYDRO_NUM_EINT_COMPONENTS .GE. 2) sumEi = sumEi + eele
  if (HYDRO_NUM_EINT_COMPONENTS .GE. 3 .AND.present(erad)) sumEi = sumEi + erad


  if (sumEi == 0.0) return

  scaleEi = eint / sumEi


  prevEi = eion
  eion = eion * scaleEi
  if (present(e1)) then
     e1 = e1 - prevEi + eion
  end if

  if (HYDRO_NUM_EINT_COMPONENTS .GE. 2) then
     prevEi = eele
     eele = eele * scaleEi
     if (present(e2)) then
        e2 = e2 - prevEi + eele
     end if
  end if

  if (HYDRO_NUM_EINT_COMPONENTS .GE. 3 .AND.present(erad)) then
     prevEi = erad
     erad = erad * scaleEi
     if (present(e3)) then
        e3 = e3 - prevEi + erad
     end if
  end if

#endif

end subroutine Hydro_recalibrateEintsForCell
