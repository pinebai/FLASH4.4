!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/gr_ptApplyBCsOneBlk
!!
!! NAME
!!  gr_ptApplyBCsOneBlk
!!
!! SYNOPSIS
!!
!!  call gr_ptApplyBCsOneBlk(integer(IN),dimension(LOW:HIGH,MDIM) :: blkLimits, &
!!                           integer(IN),dimension(LOW:HIGH,MDIM) :: blkLimitsGC, &
!!                           integer(IN)                          :: blockID)
!!
!! DESCRIPTION
!!
!! Routine to apply the external boundary conditions to a block.
!!
!! ARGUMENTS
!!  blkLimits:  Array holding the upper and lower indicies of 
!!              the interior of a block (no guardcells!)             
!!  blkLimitsGC:  Array holding the upper and lower indicies
!!                of an entire block, that is including guardcells
!!  blockID:  The block ID of the block we will apply boundary conditions to
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptApplyBCsOneBlk (blkLimits, blkLimitsGC, blockID)

#include "constants.h"
#include "Flash.h"

  use Grid_interface, ONLY : Grid_getBlkBC
  use gr_ptData, ONLY : gr_ptBuf 
  use gr_ptMapData, ONLY : gr_ptSmearLen
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: blkLimits, blkLimitsGC
  integer, intent(IN) :: blockID

  integer, dimension(LOW:HIGH, MDIM) :: faces, ignoreMe 
  integer, dimension(MDIM) :: blkSize, guard, boundaryCondition, &
       sides, guardCellID, initialCoord, newCoord
  integer :: i, j, k, axis, boundaryConditionLine
  integer :: ib, ie, jb, je, kb, ke, allcenters, iCoord, jCoord, kCoord
  logical :: coordsChange

  !No need to apply boundary conditions if we have a NGP mapping scheme.
  if (gr_ptSmearLen == 0) return


  guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1


#ifdef FLASH_GRID_UG
  call Grid_getBlkBC(blockID, faces)
#else
  call Grid_getBlkBC(blockID, ignoreMe, faces)
#endif


  allCenters = 2**NDIM
  kb=blkLimitsGC(LOW,KAXIS); ke=kb+K3D*(guard(KAXIS)-1)
  do k=LEFT_EDGE,LEFT_EDGE+K3D*(RIGHT_EDGE-1)
     jb=blkLimitsGC(LOW,JAXIS); je=jb+K2D*(guard(JAXIS)-1)
     do j=LEFT_EDGE,LEFT_EDGE+K2D*(RIGHT_EDGE-1)
        ib=blkLimitsGC(LOW,IAXIS); ie=ib+guard(IAXIS)-1
        do i=LEFT_EDGE,RIGHT_EDGE

           if ((i*j*k) /= allCenters) then !! make sure it is gCell region

              boundaryCondition(:) = NOT_BOUNDARY
              sides(:) = -1

              guardCellID(IAXIS)=i; guardCellID(JAXIS)=j; guardCellID(KAXIS)=k;

              !We have already initialised the boundary(:) array to NOT_BOUNDARY.
              !The results in this array will be changed if we detect a different
              !boundary condition.
              do axis = 1, NDIM, 1

                 if (guardCellID(axis) /= CENTER) then
                    !LEFT_EDGE corresponds to LOW side. (1+1)/2 = 1
                    !RIGHT_EDGE corresponds to HIGH side. (3+1)/2 = 2
                    sides(axis) = (guardCellID(axis) + 1) / 2
                    boundaryCondition(axis) = faces(sides(axis), axis)
                 end if

              end do


              !Perform the appropriate action depending on the boundary conditions.
              if((boundaryCondition(IAXIS) == OUTFLOW).or.&
                   (boundaryCondition(JAXIS) == OUTFLOW).or.&
                   (boundaryCondition(KAXIS) == OUTFLOW) .OR.&
                   (boundaryCondition(IAXIS) == DIODE).or.&
                   (boundaryCondition(JAXIS) == DIODE).or.&
                   (boundaryCondition(KAXIS) == DIODE) .OR.&
                   (boundaryCondition(IAXIS) == HYDROSTATIC_NVOUT).or.&
                   (boundaryCondition(JAXIS) == HYDROSTATIC_NVOUT).or.&
                   (boundaryCondition(KAXIS) == HYDROSTATIC_NVOUT) .OR.&
                   (boundaryCondition(IAXIS) == HYDROSTATIC_NVDIODE).or.&
                   (boundaryCondition(JAXIS) == HYDROSTATIC_NVDIODE).or.&
                   (boundaryCondition(KAXIS) == HYDROSTATIC_NVDIODE)) then

                 ! ---------------------------------------------------------- !
                 ! OUTFLOW/DIODE conditions need to be handled differently.   !
                 ! If we encounter an OUTFLOW or DIODE, we must zero the      !
                 ! whole guard cell section.                                                   !
                 ! ---------------------------------------------------------- !                 
                 gr_ptBuf(ib:ie, jb:je, kb:ke) = 0.0

              else

                 ! ---------------------------------------------------------- !
                 ! Otherwise we must consider REFLECTING, PERIODIC            !
                 ! (sometimes) and NOT_BOUNDARY conditions collectively.      !
                 ! i.e. we must find each guard cell's new coordinates after  !
                 ! applying *all* boundary conditions.                        !
                 ! ---------------------------------------------------------- !               
                 do kCoord = kb, ke, 1
                    do jCoord = jb, je, 1
                       do iCoord = ib, ie, 1
                          
                          initialCoord(IAXIS) = iCoord; initialCoord(JAXIS) = jCoord; initialCoord(KAXIS) = kCoord
                          coordsChange = .false.
                          

                          do axis = 1, MDIM, 1

                             if(boundaryCondition(axis) == REFLECTING .OR. &
                                  boundaryCondition(axis) == HYDROSTATIC_NVREFL) then

                                !If we have a REFLECTING boundary condition we will need 
                                !to know the line of reflection.
                                boundaryConditionLine = guard(axis) + (blkSize(axis)*(sides(axis)-1))

                                !To obtain reflected value:
                                !Find distance from line of reflection(xl) to cell position(x) = (xl - x)
                                !Add this distance to the line of reflection, then add one = (xl - x) + xl + 1
                                !This is the new coordinate value, which can be re-written as (2.xl - x + 1)
                                newCoord(axis) = (2 * boundaryConditionLine) & 
                                     - initialCoord(axis) + 1


                             else if(boundaryCondition(axis) == PERIODIC) then

                                if((faces(LOW, axis)==PERIODIC).and.(faces(HIGH, axis)==PERIODIC)) then

                                   !This can only happen if there is 1 block in the computational domain.
                                   if(sides(axis)==LOW) then
                                      newCoord(axis) = initialCoord(axis) + blkSize(axis)
                                   else if(sides(axis)==HIGH) then
                                      newCoord(axis) = initialCoord(axis) - blkSize(axis)
                                   end if

                                else

                                   !Leave the guard cell coordinates unchanged.  The guard cell will 
                                   !be communicated to another block later.
                                   newCoord(axis) = initialCoord(axis)

                                end if


                             else if(boundaryCondition(axis) == NOT_BOUNDARY) then

                                newCoord(axis) = initialCoord(axis)

                             else 

                                print *, "Boundary Condition:", boundaryCondition(axis)
                                call Driver_abortFlash("gr_ptApplyBCsOneBlk: Unrecognised boundary condition... exiting")

                             end if

                             !coordsChange becomes .true. when new and initial coordinates differ.
                             coordsChange = coordsChange .or. (newCoord(axis) /= initialCoord(axis))


                          end do  !End of loop over each axis


                          !If we find that the operation of the boundary conditions 
                          !has changed which cell recives mass accumulation then
                          !move the mass to the new cell.
                          if (coordsChange .eqv. .true.) then

                             gr_ptBuf(newCoord(IAXIS), newCoord(JAXIS), newCoord(KAXIS)) = &
                                  gr_ptBuf(newCoord(IAXIS), newCoord(JAXIS), newCoord(KAXIS)) + &
                                  gr_ptBuf(initialCoord(IAXIS), initialCoord(JAXIS), initialCoord(KAXIS))
                             
                             !Important to zero the boundary cell, otherwise its contribution 
                             !may be added twice.
                             gr_ptBuf(initialCoord(IAXIS), initialCoord(JAXIS), initialCoord(KAXIS)) = 0.0
                             
                          end if

                       end do
                    end do
                 end do


              end if  !End of if OUTFLOW or DIODE.
           end if   !End of if guard cell region.


           ib=ie+1;ie=ib+mod(i,CENTER)*blkSize(IAXIS)+mod(i+1,CENTER)*guard(IAXIS)-1
        end do
        jb=je+1;je=jb+mod(j,CENTER)*blkSize(JAXIS)+mod(j+1,CENTER)*guard(JAXIS)-1
     end do
     kb=ke+1;ke=kb+mod(k,CENTER)*blkSize(KAXIS)+mod(k+1,CENTER)*guard(KAXIS)-1
  end do

  return

end subroutine gr_ptApplyBCsOneBlk
