!!****h* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistry_interface
!!
!! This is the header file for the chemistry module that defines its
!! public interfaces.
!!***

Module PrimordialChemistry_interface
#include "constants.h"
#include "Flash.h"

  interface 
        subroutine PrimordialChemistry_computeDt (blockID,  blkLimits, blkLimitsGC, solnData, dt_check, dt_minloc)
         implicit none
         integer, intent(IN) :: blockID
         integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
         real,INTENT(INOUT)  :: dt_check
         integer,INTENT(INOUT) :: dt_minloc(5)
         real, pointer :: solnData(:,:,:,:)
        end subroutine PrimordialChemistry_computeDt
  end interface

  interface 
        subroutine PrimordialChemistry(blockCount,blockList,dt)
          integer, intent(IN) :: blockCount
          integer,dimension(blockCount), intent(IN) :: blockList
          real,intent(IN) :: dt
        end subroutine PrimordialChemistry
  end interface

  interface 
        subroutine PrimordialChemistry_finalize()
        end subroutine PrimordialChemistry_finalize
  end interface

  interface 
        subroutine PrimordialChemistry_init()
        end subroutine PrimordialChemistry_init
  end interface

end Module PrimordialChemistry_interface
