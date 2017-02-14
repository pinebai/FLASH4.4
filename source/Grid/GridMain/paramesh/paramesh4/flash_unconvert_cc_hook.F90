!!****if* source/Grid/GridMain/paramesh/paramesh4/flash_unconvert_cc_hook
!!
!! NAME
!!
!!  flash_unconvert_cc_hook
!!
!! SYNOPSIS
!!  
!!  call flash_unconvert_cc_hook(real(INOUT) :: datainout(;,:,:,:), 
!!                               integer(IN) :: nvars, 
!!                               integer(IN) :: i1,
!!                               integer(IN) :: i2,
!!                               integer(IN) :: j1,
!!                               integer(IN) :: j2,
!!                               integer(IN) :: k1,
!!                               integer(IN) :: k2, 
!!                               integer(IN) :: where, 
!!                               integer(IN) :: why, 
!!                     optional,integer(IN)  :: nlayers_in_data)
!!
!! DESCRIPTION
!!
!!  Given some block of data in DATAINOUT, loop over all of the blocks
!!  and convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity) back from corresponding conservative form
!!  (i.e., momentum) to their normal representation.
!!
!!  In other words, undo the effect of flash_convert_cc_hook.
!!
!! ARGUMENTS
!!
!!  datainout - the block data buffer
!!  nvars - number of variables, used to dimension datainout array 
!!  i1    - lower IAXIS bound 
!!  i2    - upper IAXIS bound
!!  j1    - lower JAXIS bound 
!!  j2    - upper JAXIS bound
!!  k1    - lower KAXIS bound 
!!  k2    - upper KAXIS bound
!!  where - whether to apply to interior or all data points
!!  why   - the reason why back conversion is being done; one of 
!!            gr_callReason_PROLONG   called for back conversion after data
!!                                    has been prolonged, i.e., interpolated
!!            gr_callReason_RESTRICT  called for back conversion after data
!!                                    has been restricted, i.e., averaged
!!  nlayers_in_data - the number of guardcell layers for which the datainout
!!                    array includes space (regardless of whether they
!!                    include meaningful values or not)
!!                         
!!
!! NOTES
!!
!!  The variables that are converted are the named cell-centered
!!  solution variables marked to be of type PER_MASS explicitly in a
!!  Config file.  Additionally, abundances and mass scalars are
!!  considered to be of type PER_MASS.
!!
!!  For proper functioning, DENS_VAR must not be marked as PER_MASS!
!!
!!
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_primitiveToConserve 
!!
!!
!!***

!!REORDER(4): datainout

subroutine flash_unconvert_cc_hook(datainout, nvars, i1,i2,j1,j2,k1,k2, where, why, &
     &                             nlayers_in_data)

  use Grid_data, ONLY: gr_convertToConsvdInMeshInterp, gr_vartypes, gr_anyVarToConvert
  use physicaldata, ONLY : int_gcell_on_cc
  use gr_flashHook_interfaces, hidden_Own_Name=>flash_unconvert_cc_hook

  implicit none
#include "constants.h"
#include "Flash.h"

!  real,dimension(NUNK_VARS,:,:,:),intent(INOUT) :: datainout 
  integer, intent(IN) :: nvars, i1,i2,j1,j2,k1,k2
  real,dimension(nvars, i1:i2,j1:j2,k1:k2),intent(INOUT) :: datainout
  integer, intent(IN) :: where, why
  integer, intent(IN), optional :: nlayers_in_data
  real :: dens_old_inv
  integer :: i, j, k, ivar
  integer :: il,iu,jl,ju,kl,ku, iskip, jskip, kskip
  integer :: nlayers

  if (.not. gr_convertToConsvdInMeshInterp) return
  if (.not. gr_anyVarToConvert) return

#ifdef GRID_WITH_MONOTONIC
  ! With monotonic interpolation, conversion *for prolongation* is not done
  ! here but in the interpolation routine.
  if (why .eq. gr_callReason_PROLONG) return
#endif


#ifdef DENS_VAR

  if (present(nlayers_in_data)) then
     nlayers = nlayers_in_data
  else
     nlayers = NGUARD
  end if

  if (where == gr_cells_INTERIOR) then
     iskip = nlayers
     jskip = nlayers * K2D
     kskip = nlayers * K3D
  else
     iskip = 0
     jskip = 0
     kskip = 0
  end if
  il = i1 + iskip
  iu = i2 - iskip
  jl = j1 + jskip
  ju = j2 - jskip
  kl = k1 + kskip
  ku = k2 - kskip
  if (where == gr_cells_GUARD) then
     iskip = nlayers
     jskip = nlayers * K2D
     kskip = nlayers * K3D
  end if

  do k = kl,ku
     do j = jl,ju
        do i = il,iu
           if (.not.(where == gr_cells_GUARD .and. &
                i .ge. il+iskip .and. i .le. iu-iskip .and. &
                j .ge. jl+jskip .and. j .le. ju-jskip .and. &
                k .ge. kl+kskip .and. k .le. ku-kskip)) then
              
              if (datainout(DENS_VAR,i,j,k) == 0.0) then
                 ! If unk(dens)==0, assume that for all ivars of interest --- namely,
                 ! the PER_MASS type ones --- unk(ivar)==0 held before the forward
                 ! conversion to conserved form; otherwise, the program should have
                 ! aborted already in the forward conversion. - KW
                 dens_old_inv = 0.0
              else
                 dens_old_inv = 1./datainout(DENS_VAR,i,j,k)
              end if

              do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                 if (int_gcell_on_cc(ivar)) then
                    if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
                       datainout(ivar,i,j,k) = dens_old_inv*datainout(ivar,i,j,k)
                    end if
                 end if
              end do
           end if
        end do
     end do
  end do
  
#endif

  return
end subroutine flash_unconvert_cc_hook
