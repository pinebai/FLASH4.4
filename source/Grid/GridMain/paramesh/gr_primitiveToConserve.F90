!!****if* source/Grid/GridMain/paramesh/gr_primitiveToConserve
!!
!! NAME
!!
!!  gr_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  gr_primitiveToConserve(integer(in) :: blkList(count),
!!                         integer(in) :: count)
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity) to the corresponding conservative form
!!  (i.e., momentum) if gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
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
!!  gr_conserveToPrimitive
!!
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are now conserved.  No
!!  such mechanism exists in the code yet.
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

!!REORDER(5):unk

#include "Flash.h"

subroutine gr_primitiveToConserve(blkList,count)

  use Grid_data, ONLY: gr_meshMe, gr_convertToConsvdForMeshCalls, &
                        gr_vartypes, gr_anyVarToConvert
  use Driver_interface, ONLY : Driver_abortFlash
#ifdef FLASH_GRID_PARAMESH2
#define INT_GCELL_ON_CC(IVAR) (.TRUE.)
  use physicaldata, ONLY: unk, nguard, &
       ii1 => il_bnd, ii2 => iu_bnd, &
       jj1 => jl_bnd, jj2 => ju_bnd, &
       kk1 => kl_bnd, kk2 => ku_bnd
#else
  use physicaldata, ONLY: unk, int_gcell_on_cc
  use paramesh_dimensions, ONLY: nguard, npgs, &
       ii1 => il_bnd, ii2 => iu_bnd, &
       jj1 => jl_bnd, jj2 => ju_bnd, &
       kk1 => kl_bnd, kk2 => ku_bnd

#endif
  implicit none
#include "constants.h"

  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  integer :: i, j, k,n, block, ivar
#ifdef FLASH_GRID_PARAMESH2
  integer,parameter :: npgs=1
#endif

  if (.not. gr_anyVarToConvert) return

#ifdef DENS_VAR           
  if (gr_convertToConsvdForMeshCalls) then

     do n = 1,count
        block=blkList(n)
        do k = kk1+nguard*npgs,kk2-nguard*npgs
           do j = jj1+nguard*npgs,jj2-nguard*npgs
              do i = ii1+nguard*npgs,ii2-nguard*npgs
                 if (unk(DENS_VAR,i,j,k,block) == 0.0) then
                    do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                       if (INT_GCELL_ON_CC(ivar).AND.(gr_vartypes(ivar).eq.VARTYPE_PER_MASS)) then
                          if (unk(ivar,i,j,k,block) .ne. 0.0) then
                             ! This situation would probably lead to division by zero errors in the
                             ! unk(ivar)/unk(dens) operation when converting back from conserved form later,
                             ! if we did no check. Abort if unk(ivar)!=0 and unk(dens)==0, but let
                             ! unk(ivar)==unk(dens)==0 pass. - KW
99                           format ('[gr_primitiveToConserve] PE=',I7,', ivar=',I3,', block=',I8)
                             print 99,gr_meshMe,ivar,block
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
           if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
              unk(ivar,:,:,:,block) &
                   = unk(DENS_VAR,:,:,:,block)*unk(ivar,:,:,:,block)
           end if
        enddo
     enddo

  end if
#endif

  return
end subroutine gr_primitiveToConserve
