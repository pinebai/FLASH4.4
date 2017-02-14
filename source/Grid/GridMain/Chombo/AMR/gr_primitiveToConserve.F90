!!****if* source/Grid/GridMain/Chombo/AMR/gr_primitiveToConserve
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
!!  array solnData directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

!!REORDER(4):solnData

#include "constants.h"
#include "Flash.h"

subroutine gr_primitiveToConserve(blkList,count)

  use Grid_data, ONLY: gr_meshMe, gr_convertToConsvdForMeshCalls, &
                        gr_vartypes, gr_anyVarToConvert
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: i, j, k, n, blockID, ivar

  if (.not. gr_anyVarToConvert) return

#ifdef DENS_VAR           
  if (gr_convertToConsvdForMeshCalls) then

     do n = 1,count
        blockID = blkList(n)
        call Grid_getBlkPtr(blockID, solnData)
        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 if (solnData(DENS_VAR,i,j,k) == 0.0) then
                    do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                       if (gr_vartypes(ivar).eq.VARTYPE_PER_MASS) then
                          if (solnData(ivar,i,j,k) .ne. 0.0) then
                             ! This situation would probably lead to division by zero errors in the
                             ! solnData(ivar)/solnData(dens) operation when converting back from conserved form later,
                             ! if we did no check. Abort if solnData(ivar)!=0 and solnData(dens)==0, but let
                             ! solnData(ivar)==solnData(dens)==0 pass. - KW
99                           format ('[gr_primitiveToConserve] PE=',I7,', ivar=',I3,', block=',I8)
                             print 99,gr_meshMe,ivar,blockID
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


        !This subroutine is called before Chombo fills guard cells.  This
        !means the guard cells contain invalid values.  In Chombo all
        !floating point memory is initialized to the value specied in REAL.H.

        !#  ifdef CH_USE_FLOAT /* precision */
        !#    define BASEFAB_REAL_SETVAL 1.23456789e+30
        !#  else   /* precision */
        !#    define BASEFAB_REAL_SETVAL 1.23456789e+300
        !#  endif  /* precision */

        !For FLASH (and double precision) we have the situation when 1E+300
        !is multiplied by 1E+300 and so we get a floating point overflow
        !for guard cell values.  Our solution is to convert internal
        !values from primitive to conservative.  After the guard cell fill
        !we convert all values from conservative to primitive, as now
        !guard cells also contain valid values.

        !DEV CD: Understand why 1E+300 is not overwritten with valid
        !values earlier.

        do ivar = UNK_VARS_BEGIN, UNK_VARS_END
           if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                       solnData(ivar,i,j,k) &
                            = solnData(DENS_VAR,i,j,k)*solnData(ivar,i,j,k)
                    end do
                 end do
              end do
           end if
        enddo

        call Grid_releaseBlkPtr(blockID, solnData)
     enddo

  end if
#endif

  return
end subroutine gr_primitiveToConserve
