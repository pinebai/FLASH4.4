!!****if* source/Grid/GridMain/paramesh/gr_conserveToPrimitive
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

!!REORDER(5):unk

subroutine gr_conserveToPrimitive(blkList,count,allCells)

#include "Flash.h"

  use Grid_data, ONLY : gr_smallrho,gr_smalle, gr_convertToConsvdForMeshCalls, &
                        gr_vartypes, gr_anyVarToConvert
  use physicaldata, ONLY:unk
  use tree, ONLY: nodetype
#ifdef FLASH_GRID_PARAMESH2
  use tree, ONLY: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd
#else
  use paramesh_dimensions, ONLY: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd
  use physicaldata, ONLY: no_permanent_guardcells
#endif

  implicit none
#include "constants.h"

  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  logical, intent(IN):: allCells
  real :: dens_old_inv
  integer :: i, j, k,n, block, ivar
  integer :: il,iu,jl,ju,kl,ku
  logical :: loc_allCells



  if (.not. gr_convertToConsvdForMeshCalls) return

  loc_allCells = allCells
#ifndef FLASH_GRID_PARAMESH2
  if (no_permanent_guardcells) loc_allcells = .TRUE.
#endif

  il = il_bnd
  iu = iu_bnd
  jl = jl_bnd
  ju = ju_bnd
  kl = kl_bnd
  ku = ku_bnd
  if (.not. loc_allCells) then
     il = il + NGUARD
     iu = iu - NGUARD
#if NDIM > 1
     jl = jl + NGUARD
     ju = ju - NGUARD
#endif
#if NDIM > 2
     kl = kl + NGUARD
     ku = ku - NGUARD
#endif
  end if

  do n = 1,count
     block=blkList(n)
#ifdef DENS_VAR
     do k = kl,ku
        do j = jl,ju
           do i = il,iu
              if (gr_anyVarToConvert) then
                 if (nodetype(block) .NE. ANCESTOR &
                      .OR. unk(DENS_VAR,i,j,k,block) .NE. 0.0 &
                      .OR. (iu - il) .EQ. (NXB - 1) &
                      .OR. (i .GE. il_bnd+NGUARD .AND. i .LE. iu_bnd-NGUARD &
                      .AND. j .GE. jl_bnd+NGUARD*K2D .AND. j .LE. ju_bnd-NGUARD*K2D &
                      .AND. k .GE. kl_bnd+NGUARD*K3D .AND. k .LE. ku_bnd-NGUARD*K3D)) then
                    if (unk(DENS_VAR,i,j,k,block) == 0.0) then
                       ! If unk(dens)==0, assume that for all ivars of interest --- namely,
                       ! the PER_MASS type ones --- unk(ivar)==0 held before the forward
                       ! conversion to conserved form; otherwise, the program should have
                       ! aborted already in the forward conversion. - KW
                       dens_old_inv = 0.0
                    else
                       dens_old_inv = 1./unk(DENS_VAR,i,j,k,block)
                    end if
                 else
                    dens_old_inv = 1./max(unk(DENS_VAR,i,j,k,block), gr_smallrho)
                 end if

                 ! small limits -- in case the interpolants are not monotonic
                 unk(DENS_VAR,i,j,k,block) = &
                      max(unk(DENS_VAR,i,j,k,block), gr_smallrho)

                 do ivar = UNK_VARS_BEGIN, UNK_VARS_END
                    if (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS) then
                       unk(ivar,i,j,k,block) = dens_old_inv*unk(ivar,i,j,k,block)
                    end if
                 end do
              end if

              ! small limits -- in case the interpolants are not monotonic
              unk(DENS_VAR,i,j,k,block) = &
                   max(unk(DENS_VAR,i,j,k,block), gr_smallrho)

           end do
        end do
     end do
#endif

#ifdef ENER_VAR               
     ! energy
     unk(ENER_VAR,:,:,:,block) = max(unk(ENER_VAR,:,:,:,block), gr_smalle)
#endif
#ifdef EINT_VAR
     unk(EINT_VAR,:,:,:,block) = max(unk(EINT_VAR,:,:,:,block), gr_smalle)
#endif
  end do

  return 
end subroutine gr_conserveToPrimitive
        
