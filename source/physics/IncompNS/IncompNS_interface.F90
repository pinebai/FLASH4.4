!!****h* source/physics/IncompNS/IncompNS_interface
!!
!! NAME
!!
!!  IncompNS_interface
!!
!! SYNOPSIS
!!
!!  use IncompNS_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Incompressible Navier-Stokes (INS)
!! module that defines its public interfaces.
!!
!!***

Module IncompNS_interface

  implicit none

#include "Flash.h"


  interface ! IncompNS_computeDt

    subroutine IncompNS_computeDt(ins_mindt,ins_minloc)
      implicit none
      real, intent(INOUT) :: ins_mindt
      integer, intent(INOUT) :: ins_minloc(5)
    end subroutine IncompNS_computeDt

  end interface


  interface !IncompNS
    subroutine IncompNS ( blockCount,  blockList, &
                      timeEndAdv,  dt, dtOld, &
                      sweepOrder              )
      implicit none 
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
      real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  
    end subroutine IncompNS
  end interface


  interface !IncompNS_init
    subroutine IncompNS_init (restart)
      implicit none
      logical, intent(IN) :: restart
    end subroutine IncompNS_init
  end interface


  interface !IncompNS_finalize
    subroutine IncompNS_finalize ()
      implicit none
    end subroutine IncompNS_finalize
  end interface


  interface !IncompNS_sendOutputData
    subroutine IncompNS_sendOutputData ()
      implicit none
    end subroutine IncompNS_sendOutputData
  end interface

  interface
    subroutine IncompNS_stats()
      implicit none
    end subroutine
  end interface

  interface
     subroutine IncompNS_statsIOExport(expt_flag)
       implicit none
       logical, intent(in) :: expt_flag
     end subroutine IncompNS_statsIOExport
  end interface

  interface
     subroutine IncompNS_velomg2center(blockList,blockCount)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) ::  blockCount
       integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
       !! -----------------------------------------------------
     end subroutine IncompNS_velomg2center
  end interface

end Module IncompNS_interface


