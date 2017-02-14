!!****h* source/physics/Hydro/Hydro_interface
!!
!! This is the header file for the hydro module that defines its
!! public interfaces.
!!***
Module Hydro_interface
#include "constants.h"
#include "Flash.h"

  implicit none

  interface Hydro_computeDt
     subroutine Hydro_computeDt (block_no,  &
          x, dx, uxgrid, &
          y, dy, uygrid, &
          z, dz, uzgrid, &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_check, dt_minloc, extraInfo )
       implicit none
       integer, intent(IN) :: block_no
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif
       real,INTENT(INOUT)    :: dt_check
       integer,INTENT(INOUT)    :: dt_minloc(5)
       real, pointer :: solnData(:,:,:,:) 
       real, OPTIONAL,intent(INOUT) :: extraInfo
     end subroutine Hydro_computeDt
  end interface

  interface
     subroutine Hydro_consolidateCFL
     end subroutine Hydro_consolidateCFL
  end interface

  interface Hydro
     subroutine Hydro ( blockCount, blockList, &
                        timeEndAdv, dt, dtOld,  &
                        sweepOrder )

       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(blockcount) :: blockList
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
     end subroutine Hydro
  end interface

  interface Hydro_init
     subroutine Hydro_init()
     end subroutine Hydro_init
  end interface

  interface Hydro_finalize
     subroutine Hydro_finalize()
     end subroutine Hydro_finalize
  end interface

  interface Hydro_detectShock
     subroutine Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, &
          guardCells, &
          primaryCoord,secondCoord,thirdCoord)
       
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
       integer, intent(IN) :: guardCells(MDIM)
       real, pointer,dimension(:,:,:,:) :: solnData
#ifdef FIXEDBLOCKSIZE
       real, intent(out),dimension(GRID_ILO_GC:GRID_IHI_GC,&
                                   GRID_JLO_GC:GRID_JHI_GC,&
                                   GRID_KLO_GC:GRID_KHI_GC):: shock
       real,intent(IN),dimension(GRID_ILO_GC:GRID_IHI_GC) :: primaryCoord
       real,intent(IN),dimension(GRID_JLO_GC:GRID_JHI_GC) :: secondCoord
       real,intent(IN),dimension(GRID_KLO_GC:GRID_KHI_GC) :: thirdCoord
#else
       real,intent(out), dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
       real,intent(IN),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: primaryCoord
       real,intent(IN),dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: secondCoord
       real,intent(IN),dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: thirdCoord
#endif

     end subroutine Hydro_detectShock
  end interface
  
  interface Hydro_shockStrength
     subroutine Hydro_shockStrength(solnData, shock, blkLimits, blkLimitsGC, &
          guardCells, &
          primaryCoord,secondCoord,thirdCoord, &
          threshold, mode)
          implicit none
       integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
       integer, intent(IN) :: guardCells(MDIM)
       real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
       real, intent(inout),dimension(GRID_ILO_GC:GRID_IHI_GC,&
            GRID_JLO_GC:GRID_JHI_GC,&
            GRID_KLO_GC:GRID_KHI_GC):: shock
       real,intent(IN),dimension(GRID_ILO_GC:GRID_IHI_GC) :: primaryCoord
       real,intent(IN),dimension(GRID_JLO_GC:GRID_JHI_GC) :: secondCoord
       real,intent(IN),dimension(GRID_KLO_GC:GRID_KHI_GC) :: thirdCoord
#else
       real,intent(inout),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
       real,intent(IN),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: primaryCoord
       real,intent(IN),dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: secondCoord
       real,intent(IN),dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: thirdCoord
#endif
       real, intent(IN) :: threshold
       integer, intent(IN) :: mode
     end subroutine Hydro_shockStrength
  end interface


  interface Hydro_sendOutputData
     subroutine Hydro_sendOutputData()
       
     end subroutine Hydro_sendOutputData
  end interface

  interface Hydro_recalibrateEints
     subroutine Hydro_recalibrateEints(range,blockID)
       implicit none
       integer, dimension(2,MDIM), intent(in) :: range
       integer,intent(in) :: blockID
     end subroutine Hydro_recalibrateEints
     subroutine Hydro_recalibrateEintsForCell(eint,eion,eele,erad,e1,e2,e3)
       implicit none
       real,intent(in)    :: eint
       real,intent(INOUT) :: eion,eele
       real,intent(INOUT),OPTIONAL :: erad
       real,intent(INOUT),OPTIONAL :: e1,e2,e3
     end subroutine Hydro_recalibrateEintsForCell
  end interface

  interface
     logical function Hydro_gravPotIsAlreadyUpdated()
       implicit none
     end function Hydro_gravPotIsAlreadyUpdated
  end interface

  interface
     subroutine Hydro_mapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
          axis,face,idest)
       implicit none
       integer, intent(OUT) :: bcTypeToApply
       integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
       integer,intent(IN),OPTIONAL:: idest
     end subroutine Hydro_mapBcType
  end interface

end Module Hydro_interface
