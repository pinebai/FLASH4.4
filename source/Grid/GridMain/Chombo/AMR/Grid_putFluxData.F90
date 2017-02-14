!!****if* source/Grid/GridMain/Chombo/AMR/Grid_putFluxData
!!
!! NAME
!!  Grid_putFluxData
!!
!! SYNOPSIS
!!
!!
!!  Grid_putFluxData(integer(IN) :: blockID, 
!!                   integer(IN) :: axis,
!!                   real(IN)    :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots,
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!
!!  Put the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used with adaptive mesh
!!  since fluxs calculated by the two blocks that are at fine/coarse boundary have 
!!  different accuracy. The fluxes calculated by individual blocks are reported to 
!!  the Grid through this call. Once that is done, a call to Grid_conserveFluxes 
!!  applies the flux conservation algorithm to make it consistent across the fine/coarse 
!!  boundaries.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to put fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.
!!
!!
!!  dataSize : integer array specifying the dimensions for the array, fluxes
!!
!!             dataSize (1) holds the number of cells provided in the i direction
!!
!!             dataSize (2) holds the number of cells provided in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells provided in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells, and the  fluxes must be correct for 
!!             the interior cells of the block, as this interface does not know which 
!!             cell fluxes the Grid will need to store.
!!
!!  pressureSlots: If present and greater than zero, this indicates one flux variable
!!                 in the fluxes array that may need special handling because it
!!                 really scales like a flux; normally this would be pressure,
!!                 but it could be another flux variable that the caller keeps in
!!                 flux density form. Ignored in this implementation since with
!!                 Paramesh2, all "fluxes" are assumed to be give as flux densities
!!                 anyway as far as the Grid unit is concerned.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlot argument.
!!                 Ignored in this implementation since with
!!                 Paramesh2, all "fluxes" are assumed to be give as flux densities
!!                 anyway as far as the Grid unit is concerned.
!!
!! NOTES 
!!
!!   This implementation is specific to Chombo.
!!
!!***

!!REORDER(4): fluxes, c_fluxes

#include "constants.h"
#include "Flash.h"

subroutine Grid_putFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkRefineLevel
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getDt
  use iso_c_binding, ONLY : c_loc, c_double, c_int
  use chombo_f_c_interface, ONLY : ch_put_flux_data
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Grid_data, ONLY : gr_verbosity

  implicit none
  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(IN), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL :: pressureSlots
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

  real, target, allocatable, dimension(:,:,:,:) :: c_fluxes
  real :: dt
  real(c_double) :: c_dt
  integer(c_int) :: c_blockID, c_axis
  integer, dimension(NUNK_VARS) :: solnFluxMap
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: fSize, skip
  integer, parameter :: no_map = -1  !Must be negative to avoid clashes
  integer :: uIdx, fIdx
  integer :: logUnit, i, j, k, v, refineLevel
  logical :: debugFluxes

#ifdef DEBUG_FLUXES
  debugFluxes = .true.
#else
  debugFluxes = .false.
#endif


#if NFLUXES > 0

  if (NPROP_FLUX /= 7) then
     !We should extend the setup script so that the mapping
     !between solution and flux properties is dynamic.
     call Driver_abortFlash("We have hard-coded the flux to solution mapping")
  else

     call Driver_getDt(dt)
     c_dt = dt
     c_blockID = blockID
     c_axis = axis


     solnFluxMap(:) = no_map


     !Testing for the existance of property fluxes is a hack so that
     !EOS unit test still works.
#if NPROP_FLUX > 0
#if NPROP_FLUX == 7 && defined(FLASH_HYDRO_PPM)
     solnFluxMap(DENS_VAR) = RHO_FLUX
     solnFluxMap(EINT_VAR) = EINT_FLUX
     solnFluxMap(ENER_VAR) = E_FLUX
     solnFluxMap(PRES_VAR) = P_FLUX

     if (axis == IAXIS) then
        solnFluxMap(VELX_VAR) = U_FLUX
        solnFluxMap(VELY_VAR) = UT_FLUX
        solnFluxMap(VELZ_VAR) = UTT_FLUX
     else if (axis == JAXIS) then
        solnFluxMap(VELX_VAR) = UT_FLUX
        solnFluxMap(VELY_VAR) = U_FLUX
        solnFluxMap(VELZ_VAR) = UTT_FLUX
     else if (axis == KAXIS) then
        solnFluxMap(VELX_VAR) = UT_FLUX
        solnFluxMap(VELY_VAR) = UTT_FLUX
        solnFluxMap(VELZ_VAR) = U_FLUX
     else
        call Driver_abortFlash("Unknown axis")
     end if
#else
     !We should extend the setup script so that the mapping
     !between solution and flux properties is dynamic.
     call Driver_abortFlash("We have hard-coded the flux to solution mapping "//&
          "for split hydro only")
#endif
#endif


     if (NUNK_VARS > NPROP_VARS) then
        !We have species or mass scalars to consider.  The mapping
        !betweeen solution and fluxes is simple for species and mass scalars.
        fIdx = PROP_FLUX_END + 1
        do uIdx = PROP_VARS_END+1, NUNK_VARS
           solnFluxMap(uIdx) = fIdx
           fIdx = fIdx + 1
        end do
     end if


     !If the cell-centered block has 8 internal cells and 2*4 guard cells then
     !the internal cells are at indices 5 to 12.  For face-centered fluxes
     !we need to add 1 to the internal cells along the specified axis.  This
     !means that flux data will exist at internal cells 5 to 13 for one axis.
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     blkLimits(HIGH,axis) = blkLimits(HIGH,axis) + 1


     !We can optimize the frequent allocation/deallocation later.
     fSize(1:MDIM) = 1
     fSize(1:NDIM) = blkLimits(HIGH,1:NDIM) - blkLimits(LOW,1:NDIM) + 1
     allocate(c_fluxes(NUNK_VARS,fSize(IAXIS),fSize(JAXIS),fSize(KAXIS)))

     do uIdx = 1, NUNK_VARS
        if (no_map == solnFluxMap(uIdx)) then
           c_fluxes(uIdx, &
                1:fSize(IAXIS), &
                1:fSize(JAXIS), &
                1:fSize(KAXIS)) = 0.0
        else
           fIdx = solnFluxMap(uIdx)
           c_fluxes(uIdx, &
                1:fSize(IAXIS), &
                1:fSize(JAXIS), &
                1:fSize(KAXIS)) = fluxes(fIdx, &
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
        end if
     end do

     c_fluxes(PRES_VAR,:,:,:) = 0.0
     c_fluxes(EINT_VAR,:,:,:) = 0.0


     if (gr_verbosity >= 4) then
        !Presenting data this way allows us to compare against the output
        !from poutCoarseRegisters() and poutFineRegisters().
        call Logfile_open(logUnit,.true.)
        skip = 1; skip(axis) = fSize(axis)-1
        call Grid_getBlkRefineLevel(blockID, refineLevel)
        write (logUnit,"(a,i3,a,i3,a,i4,a)") "Fluxes for axis", axis-1, &
             ", level", refineLevel-1, ", block", blockID-1, &
             " (each column is a cell-centered variable in Flash.h order)"

        do k = 1, fSize(KAXIS), skip(KAXIS)
           do j = 1, fSize(JAXIS), skip(JAXIS)
              do i = 1, fSize(IAXIS), skip(IAXIS)
                 if (NDIM == 1) then
                    write (logUnit,"(a,i2,a)",advance="no") &
                         "(", i-1, ")"
                 else if (NDIM == 2) then
                    write (logUnit,"(a,i2,a,i2,a)",advance="no") &
                         "(", i-1, ",", j-1, ")"
                 else if (NDIM == 3) then
                    write (logUnit,"(a,i2,a,i2,a,i2,a)",advance="no") &
                         "(", i-1, ",", j-1, ",", k-1, ")"
                 else
                    call Driver_abortFlash("Invalid dimensionality")
                 end if
                    
                 do v = 1, NUNK_VARS
                    write (logUnit,"(ES12.3)",advance="no") c_fluxes(v,i,j,k)
                 end do
                 write (logUnit,*)
              end do
           end do
        end do
        write (logUnit,*)
        call Logfile_close(.true.)
     end if


     if (debugFluxes) then
        !Compare with pout.0 output.  Need to change the scientfic notation
        !character E to e and remove the end space character from each line.
        !sed 's/E/e/g;s/ $//' 2blast_chombo_amr_1d0.log > compare_to_pout.0
        !diff pout.0 compare_to_pout.0
        call Logfile_open(logUnit,.true.)
        do v = 1, NUNK_VARS
           do k = 1, fSize(KAXIS)
              do j = 1, fSize(JAXIS)
                 do i = 1, fSize(IAXIS)
                    write (logUnit,"(ES12.3)",advance="no") c_fluxes(v,i,j,k)
                 end do
              end do
           end do
           write (logUnit,*)
        end do
        write (logUnit,*)
        call Logfile_close(.true.)
     end if

     !Pass a pointer to the first element of the c_fluxes array.  The entire
     !extent of allocated data will initialize the Chombo FArrayBox for fluxes
     !directly (there will be no memory skips).
     call ch_put_flux_data(c_loc(c_fluxes(1,1,1,1)), c_dt, c_blockID, c_axis)
     deallocate(c_fluxes)

  end if
#endif

end subroutine Grid_putFluxData
