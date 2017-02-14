!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftGetProcGrid
!!
!! NAME
!!
!! gr_pfftGetProcGrid
!!
!! SYNOPSIS
!!  
!! gr_pfftGetProcGrid(integer(IN) :: dims, &
!!                    integer(IN) :: pfftMyPE, &
!!                    integer(IN) :: pfftNumProcs, &
!!                    integer(IN) :: pfftGlobalLen(MDIM), &
!!                    integer(OUT) :: pfftProcGrid(MDIM))
!!
!! DESCRIPTION
!!
!! Gets the pencil processor grid.
!!
!! ARGUMENTS
!!
!! dims - The dimensionality of the simulation.
!! pfftMyPE - My PE.
!! pfftNumProcs - The number of processors we must use.
!! pfftGlobalLen - An array containing the number of grid points
!!                 in the global pencil grid.
!! pfftProcGrid - An array which will contain the processor grid.
!!
!! NOTES
!!
!! Each processor will reach the same conclusion, i.e. the same
!! calculations are performed on each processor.  The subroutine 
!! is relatively quick to execute, so this strategy is appropriate.
!!
!!***
subroutine gr_pfftGetProcGrid(dims, pfftMyPE, pfftNumProcs, &
     pfftGlobalLen, pfftProcGrid)
#include "constants.h"
  use gr_pfftInterfaceTypeDecl, ONLY : PossibleGrid_t
  use gr_pfftInterface, ONLY : gr_pfftMakePencilIn3dSpace, &
       gr_pfftFnArgHardConstraint, gr_pfftFnArgMediumConstraint, &
       gr_pfftFnArgEasyConstraint
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshNumProcs

  implicit none
  integer, intent(IN) :: dims, pfftMyPE, pfftNumProcs
  integer, dimension(MDIM), intent(IN) :: pfftGlobalLen
  integer, dimension(MDIM), intent(OUT) :: pfftProcGrid
  type(PossibleGrid_t) :: pencilGrid


  pfftProcGrid(1:MDIM) = 1
  if ((pfftNumProcs < 1) .or. (pfftNumProcs > gr_meshNumProcs)) then
     call Driver_abortFlash("[gr_pfftGetProcGrid]: Number of procs invalid!")
  end if

  select case (dims)
  case (1)
     if (pfftNumProcs /= 1) then
        call Driver_abortFlash &
             ("Pfft needs to have one dimension in processor")
     end if

  case (2)
     pfftProcGrid(JAXIS) = pfftNumProcs

  case (3)
     !We attempt to generate the PFFT grid according to various levels of 
     !constraints. The global grid size is a multiple of the PFFT processor 
     !grid when we satisfy the hard constraint.  This is a nice scenario to 
     !work with.  Things become more difficult when different amounts of 
     !data move to different processors, as is the case with the medium and 
     !easy constraints.

     pencilGrid = gr_pfftMakePencilIn3dSpace &
          (pfftGlobalLen, pfftNumProcs, gr_pfftFnArgHardConstraint)
     if ( (pencilGrid % jProcs == NONEXISTENT) .or. &
          (pencilGrid % kProcs == NONEXISTENT) ) then

        !Warn the user that this is a difficult grid to work with.
        pencilGrid = gr_pfftMakePencilIn3dSpace &
             (pfftGlobalLen, pfftNumProcs, gr_pfftFnArgMediumConstraint)
        if ( (pencilGrid % jProcs == NONEXISTENT) .or. &
             (pencilGrid % kProcs == NONEXISTENT) ) then

           !Failure to satisfy easy constraint should never happen!
           pencilGrid = gr_pfftMakePencilIn3dSpace &
                (pfftGlobalLen, pfftNumProcs, gr_pfftFnArgEasyConstraint)
           if ( (pencilGrid % jProcs == NONEXISTENT) .or. &
                (pencilGrid % kProcs == NONEXISTENT) ) then
              call Driver_abortFlash &
                   ("[gr_pfftGetProcGrid]: Can't generate a valid PFFT grid???")
           else
              
              if (pfftMyPE == 0) then
                 print *, "[gr_pfftGetProcGrid]: WARNING. " // &
                      "We have a VERY difficult domain for PFFT!"
              end if
           end if

        else
           if (pfftMyPE == 0) then
              print *, "[gr_pfftGetProcGrid]: WARNING. " // &
                   "We have a difficult domain for PFFT!"
           end if
        end if
     end if

     pfftProcGrid(JAXIS) = pencilGrid % jProcs
     pfftProcGrid(KAXIS) = pencilGrid % kProcs

  case default
     call Driver_abortFlash("Invalid dimensionality")

  end select

end subroutine gr_pfftGetProcGrid
