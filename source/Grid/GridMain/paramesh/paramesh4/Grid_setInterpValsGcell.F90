!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_setInterpValsGcell
!!
!!
!! NAME
!!
!!  Grid_setInterpValsGcell
!!
!!
!! SYNOPSIS
!!
!!  call Grid_setInterpValsGcell(logical(IN) :: setval)
!!
!!
!! DESCRIPTION
!!
!! Sets interpolation values for guardcell filling within INS fractional
!! step method.
!!
!!
!! ARGUMENTS
!!
!! setval = .true.  set values 
!!          .false. restore all interpolation-restriction values to quadratic
!!
!! NOTES
!!
!!  The current functionality is specifically for use by the IncompNS unit.
!!
!! DEV: This should be generalized, to be useful for other code than IncompNS!
!!***

subroutine Grid_setInterpValsGcell(setval)

#include "Flash.h"

  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work                           
#endif    

  implicit none

  logical, intent(IN) :: setval

  integer, parameter :: intval      = 2 ! Quadratic interpolation and restriction
  integer, parameter :: intfaceval = 32 ! Quadratic restriction + linear restriction
                                        ! in block face.

!  ------------------------------------------------------------------------

#ifdef FLASH_GRID_PARAMESH

  if (setval) then ! Set Gcell interpolation and restriction values:

     ! Work interpolation:
     interp_mask_work(:)= intval;

     ! Cell centered valiables interpolation and Restriction:
     interp_mask_unk(:)     = intval;   
     interp_mask_unk_res(:) = intval;
  
     ! Face Centered variables interpolation and Restriction:
     interp_mask_facex(:) = intval; interp_mask_facex_res(:) = intval;
     interp_mask_facey(:) = intval; interp_mask_facey_res(:) = intval;
     interp_mask_facez(:) = intval; interp_mask_facez_res(:) = intval;

     ! For velocities linear restriction on the face:
#ifdef VELC_FACE_VAR
     interp_mask_facex_res(VELC_FACE_VAR) = intfaceval
     interp_mask_facey_res(VELC_FACE_VAR) = intfaceval
     interp_mask_facez_res(VELC_FACE_VAR) = intfaceval
#endif

  else ! Restore velocity restriction values to quadratic:

     ! For velocities linear restriction on the face:
#ifdef VELC_FACE_VAR
     interp_mask_facex_res(VELC_FACE_VAR) = intval
     interp_mask_facey_res(VELC_FACE_VAR) = intval
     interp_mask_facez_res(VELC_FACE_VAR) = intval
#endif

  endif
  
#endif

  return

end subroutine Grid_setInterpValsGcell

