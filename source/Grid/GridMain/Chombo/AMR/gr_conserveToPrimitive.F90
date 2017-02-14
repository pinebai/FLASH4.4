!!****if* source/Grid/GridMain/Chombo/AMR/gr_conserveToPrimitive
!!
!! NAME
!!
!!  gr_conserveToPrimitive
!!
!!
!! SYNOPSIS
!!
!!  gr_conserveToPrimitive(integer(in) :: blkList(count),
!!                         integer(in) :: count,
!!                         logical(in) :: allCells)
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity) from the corresponding conservative form
!!  (i.e., momentum) back to the normal PER_MASS form  if 
!!  gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!  Additionally,
!!   - energies (ENER_VAR and EINT_VAR) are forced to be .ge. gr_smalle,
!!   - the density (DENS_VAR) is forced to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
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
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are no longer conserved.  No
!!  such mechanism exists in the code yet.
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

!!REORDER(4):solnData

#include "constants.h"
#include "Flash.h"

subroutine gr_conserveToPrimitive(blkList,count,allCells)

  use Grid_data, ONLY : gr_smallrho,gr_smalle, gr_convertToConsvdForMeshCalls, &
                        gr_vartypes, gr_anyVarToConvert
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits

  implicit none
  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  logical, intent(IN):: allCells
  real :: dens_old_inv
  integer :: i, j, k,n, blockID, ivar
  integer :: il,iu,jl,ju,kl,ku
  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC, lim

  if (.not. gr_convertToConsvdForMeshCalls) return


  do n = 1,count

     blockID = blkList(n)
     call Grid_getBlkPtr(blockID, solnData)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     if (allCells) then
        ! Get the entire block region.
        lim = blkLimitsGC
     else
        ! Get the internal block region.
        lim = blkLimits
     end if

     il = lim(LOW,IAXIS)
     iu = lim(HIGH,IAXIS)
     jl = lim(LOW,JAXIS)
     ju = lim(HIGH,JAXIS)
     kl = lim(LOW,KAXIS)
     ku = lim(HIGH,KAXIS)


#ifdef DENS_VAR
     do k = kl,ku
        do j = jl,ju
           do i = il,iu
              if (gr_anyVarToConvert) then
                 !CD: The Paramesh version of this subroutine is much more complicated
                 !because it includes logic to avoid dens = 0.0 in ANCESTOR blocks.
                 if (solnData(DENS_VAR,i,j,k) == 0.0) then
                    ! If solnData(dens)==0, assume that for all ivars of interest --- namely,
                    ! the PER_MASS type ones --- solnData(ivar)==0 held before the forward
                    ! conversion to conserved form; otherwise, the program should have
                    ! aborted already in the forward conversion. - KW
                    dens_old_inv = 0.0
                 else
                    dens_old_inv = 1./solnData(DENS_VAR,i,j,k)
                 end if

                 ! small limits -- in case the interpolants are not monotonic
                 solnData(DENS_VAR,i,j,k) = &
                      max(solnData(DENS_VAR,i,j,k), gr_smallrho)

                 do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                    if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
                       solnData(ivar,i,j,k) = dens_old_inv*solnData(ivar,i,j,k)
                    end if
                 end do
              end if

              ! small limits -- in case the interpolants are not monotonic
              solnData(DENS_VAR,i,j,k) = &
                   max(solnData(DENS_VAR,i,j,k), gr_smallrho)

           end do
        end do
     end do
#endif

#ifdef ENER_VAR
     ! energy
     solnData(ENER_VAR,:,:,:) = max(solnData(ENER_VAR,:,:,:), gr_smalle)
#endif
#ifdef EINT_VAR
     solnData(EINT_VAR,:,:,:) = max(solnData(EINT_VAR,:,:,:), gr_smalle)
#endif

     call Grid_releaseBlkPtr(blockID, solnData)
  end do

  return
end subroutine gr_conserveToPrimitive
