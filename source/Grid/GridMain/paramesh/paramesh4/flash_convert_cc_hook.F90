!!****if* source/Grid/GridMain/paramesh/paramesh4/flash_convert_cc_hook
!!
!! NAME
!!
!!  flash_convert_cc_hook
!!
!!
!! SYNOPSIS
!!
!!  call flash_convert_cc_hook(real(INOUT) :: datainout(:,:,:,:),
!!                             integer(IN) :: nvars, 
!!                             integer(IN) :: i1,
!!                             integer(IN) :: i2,
!!                             integer(IN) :: j1,
!!                             integer(IN) :: j2,
!!                             integer(IN) :: k1,
!!                             integer(IN) :: k2, 
!!                             integer(IN) :: why)
!!
!! DESCRIPTION
!!
!!  Given some block of data in DATAINOUT (usually copied or extracted
!!  from the global solution array UNK), loop over all variables and
!!  convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity, specific energies, mass fractions) to the
!!  corresponding conservative form (i.e., momentum, energy density,
!!  partial mass densities, respectively).
!!
!!
!! ARGUMENTS
!! 
!!  datainout - one block's worth of unk1 data
!!  nvars - number of variables, used to dimension datainout array 
!!  i1    - lower IAXIS bound 
!!  i2    - upper IAXIS bound
!!  j1    - lower JAXIS bound 
!!  j2    - upper JAXIS bound
!!  k1    - lower KAXIS bound 
!!  k2    - upper KAXIS bound
!!  why   - the reason why conversion is being done; one of 
!!            gr_callReason_PROLONG   called for conversion before data
!!                                    is prolonged, i.e., interpolated
!!            gr_callReason_RESTRICT  called for conversion before data
!!                                    is restricted, i.e., summarized
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
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_primitiveToConserve
!!
!!***

!!REORDER(4):datainout

subroutine flash_convert_cc_hook(datainout, nvars, i1,i2,j1,j2,k1,k2, why)

  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_data, ONLY: gr_convertToConsvdInMeshInterp, gr_vartypes, gr_anyVarToConvert, &
       gr_meshMe
  use physicaldata, ONLY : int_gcell_on_cc
  use gr_flashHook_interfaces, hidden_Own_Name=>flash_convert_cc_hook

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: nvars, i1,i2,j1,j2,k1,k2
  real,dimension(nvars, i1:i2,j1:j2,k1:k2),intent(INOUT) :: datainout
  integer, intent(IN) :: why

  integer :: i, j, k, ivar
  integer :: ii1,ii2,jj1,jj2,kk1,kk2
  if (.not. gr_convertToConsvdInMeshInterp) return
  if (.not. gr_anyVarToConvert) return

#ifdef GRID_WITH_MONOTONIC
  ! With monotonic interpolation, conversion *for prolongation* is not done
  ! here but in the interpolation routine.
  if (why .eq. gr_callReason_PROLONG) return
#endif


#ifdef DENS_VAR           
  if (why .eq. gr_callReason_RESTRICT) then
     ii1 = i1+NGUARD
     ii2 = i2-NGUARD
     jj1 = j1+K2D*NGUARD
     jj2 = j2-K2D*NGUARD
     kk1 = k1+K3D*NGUARD
     kk2 = k2-K3D*NGUARD
  else 
     ii1 = i1
     ii2 = i2
     jj1 = j1
     jj2 = j2
     kk1 = k1
     kk2 = k2
  end if

  do k = kk1,kk2
     do j = jj1,jj2
        do i = ii1,ii2
           if (datainout(DENS_VAR,i,j,k) == 0.0) then
              do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                 if (int_gcell_on_cc(ivar).AND.(gr_vartypes(ivar).eq.VARTYPE_PER_MASS)) then
                    if (datainout(ivar,i,j,k) .ne. 0.0) then
                       ! This situation would probably lead to division by zero errors in the
                       ! unk(ivar)/unk(dens) operation when converting back from conserved form later,
                       ! if we did no check. Abort if unk(ivar)!=0 and unk(dens)==0, but let
                       ! unk(ivar)==unk(dens)==0 pass. - KW
99                     format ('[flash_convert_cc_hook] PE=',I7,', ivar=',I3,', why=',I1,', i,j,k=',3I5)
                       print 99,gr_meshMe,ivar,why,i,j,k
                       print*,'Trying to convert non-zero mass-specific variable to per-volume form, but dens is zero!'
                       call Driver_abortFlash &
                            ('Trying to convert non-zero mass-specific variable to per-volume form, but dens is zero!')
                    end if
                 end if
              end do
           end if
        end do
     end do
  end do

  do ivar = UNK_VARS_BEGIN, UNK_VARS_END
     if (int_gcell_on_cc(ivar)) then
        if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
           datainout(ivar,ii1:ii2,jj1:jj2,kk1:kk2) &
                = datainout(DENS_VAR,ii1:ii2,jj1:jj2,kk1:kk2)*datainout(ivar,ii1:ii2,jj1:jj2,kk1:kk2)
        end if
     end if
  enddo
#endif

  return
end subroutine flash_convert_cc_hook
