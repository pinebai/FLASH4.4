!!****if* source/Grid/GridMain/paramesh/gr_setGcFillNLayers
!!
!! NAME
!!
!!  gr_setGcFillNLayers
!!
!! SYNOPSIS
!!
!!  gr_setGcFillNLayers(integer(OUT)  :: layers(MDIM),
!!                      integer(IN)  :: idir,
!!                      integer(IN)  :: guard,
!!                      integer(IN),OPTIONAL  :: minlayers,
!!                      integer(OUT),OPTIONAL :: returnLayers(MDIM))
!!
!! DESCRIPTION
!!
!!   Given a caller's required number of guard cell layers that need
!!   to be filled, this routine determines the values of nlayers{x,y,z}
!!   in the call to amr_guardcell.
!!
!!   The numbers of layers requested from Paramesh are alway at
!!   least as many as those required by the caller, but may be
!!   larger for two reasons:
!!    o  to ensure that parent blocks have enough guard cells with
!!       valid data to allow interpolation into their chiildren's
!!       guard cells (this depends on the interpolation method's 
!!       stencil); and/or
!!    o  to ensure that alignment reqirements of the interpolation
!!       method are satisfied ("monotonic" interpolation requires
!!       the child cell ranges to start at odd indices).
!!
!!   This subroutine is meant to be called from Grid_fillGuardCells.
!!
!! ARGUMENTS
!!
!!   layers : returns value to pass to amr_guardcell as nlayersx/y/z argument
!!
!!
!!   idir :     selects a direction (or not) - interpreted as for Grid_fillGuardCells.
!!
!!   guard :    max number of layers possible - should probably always be the same
!!              as NGUARD.
!!
!!   minlayers : selects (or not) minimum number of guard cells *that the caller needs*,
!!               interpreted as for Grid_fillGuardCells.
!!
!!   returnLayers: if present, this array will on return simply contain the
!!               requested number of layers for each direction.
!!
!! NOTES
!!
!!   Currently this subroutine is only used in Paramesh 3 ff. versions of
!!   Grid_fillGuardCells and in Multigrid Poisson solvers.
!!
!! SEE ALSO
!!
!!   Grid_fillGuardCells
!!
!!***

subroutine gr_setGcFillNLayers(layers, idir, guard, minLayers, returnLayers)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_intpolStencilWidth

#include "Flash.h"
#include "constants.h"

  implicit none
  integer,dimension(MDIM), intent(OUT) :: layers
  integer, intent(IN)  :: idir, guard
  integer, intent(IN),OPTIONAL  :: minLayers
  integer,intent(OUT),OPTIONAL  :: returnLayers(MDIM)

  integer :: nlayers_transverse, nlayers_transverse_parents, retnlayers_transverse


  if(idir==ALLDIR .and. .not. present(minLayers))then
     layers=guard
     returnLayers=guard
  else

     ! This should work for all combinations of idir and minLayers.
     ! Basically, the minLayers requested by the caller is upped to a "safe"
     ! value so that enough layers of guard cells are exchanged to ensure
     ! that where interpolation takes place, parent blocks have enough valid
     ! cells.
     nlayers_transverse = 0
     if (present(minLayers)) then
        if (minLayers>0) nlayers_transverse = minLayers
     end if
     nlayers_transverse_parents = ((nlayers_transverse + 1) / 2) + gr_intpolStencilWidth
     nlayers_transverse = max(nlayers_transverse,nlayers_transverse_parents)
#ifdef GRID_GC_LAYERS_ALWAYS_EVEN
     nlayers_transverse = ((nlayers_transverse + 1) / 2)  * 2
#endif

     select case (idir)
     case(IAXIS)
        layers(IAXIS)=guard
        layers(JAXIS)=nlayers_transverse
        layers(KAXIS)=nlayers_transverse
     case(JAXIS)
        layers(IAXIS)=nlayers_transverse
        layers(JAXIS)=guard
        layers(KAXIS)=nlayers_transverse
     case(KAXIS)
        layers(IAXIS)=nlayers_transverse
        layers(JAXIS)=nlayers_transverse
        layers(KAXIS)=guard
     case(ALLDIR)
        layers(IAXIS)=nlayers_transverse
        layers(JAXIS)=nlayers_transverse
        layers(KAXIS)=nlayers_transverse
     case default
        print*,"[gr_setGcFillNLayers] Wrong direction specification in Grid_fillGuardCells: idir=",idir
        call Driver_abortFlash("[gr_setGcFillNLayers] Wrong direction specification in Grid_fillGuardCells")
     end select

     if (present(returnLayers)) then
        retnlayers_transverse = 0
        if (present(minLayers)) then
           if (minLayers>0) retnlayers_transverse = minLayers
        end if
        select case (idir)
        case(IAXIS)
           returnLayers(IAXIS)=guard
           returnLayers(JAXIS)=retnlayers_transverse
           returnLayers(KAXIS)=retnlayers_transverse
        case(JAXIS)
           returnLayers(IAXIS)=retnlayers_transverse
           returnLayers(JAXIS)=guard
           returnLayers(KAXIS)=retnlayers_transverse
        case(KAXIS)
           returnLayers(IAXIS)=retnlayers_transverse
           returnLayers(JAXIS)=retnlayers_transverse
           returnLayers(KAXIS)=guard
        case(ALLDIR)
           returnLayers(IAXIS)=retnlayers_transverse
           returnLayers(JAXIS)=retnlayers_transverse
           returnLayers(KAXIS)=retnlayers_transverse
        end select
     end if
  end if


end subroutine gr_setGcFillNLayers
