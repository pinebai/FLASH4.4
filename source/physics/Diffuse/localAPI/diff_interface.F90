!!****ih* source/physics/Diffuse/localAPI/diff_interface
!!
!! NAME
!!   diff_interface
!!
!! SYNOPSIS
!!   use diff_interface,ONLY: diff_advanceTherm
!!
!! DESCRIPTION
!! This is the header file for the diffuse module that defines some of its
!! internal routine interfaces.
!!***
Module diff_interface
#include "constants.h"
#include "Flash.h"


  interface
     subroutine diff_advanceTherm(blockCount,blockList,dt,pass)
       implicit none
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       real,intent(in) :: dt
       integer, OPTIONAL, intent(IN):: pass
     end subroutine diff_advanceTherm
  end interface

  interface
     subroutine diff_computeAX (blockID, blkLimits, blkLimitsGC, iVar, iFactorA, dt, theta, AX,iFactorC, iFactorD)
       implicit none
       integer, intent(IN) :: blockID
       integer, dimension(2,MDIM),intent(IN) :: blkLimits
       integer, dimension(2,MDIM),intent(IN) :: blkLimitsGC
       integer,intent(IN) :: iVar
       integer,intent(IN) :: iFactorA
       real,intent(IN):: dt
       real,intent(IN):: theta   
       real,intent(OUT) :: AX (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       integer,intent(IN),OPTIONAL :: iFactorC
       integer,intent(IN),OPTIONAL :: iFactorD
     end subroutine diff_computeAX
  end interface

  interface
     subroutine diff_computeblkAMat(blockID, NZ, N, dt, theta, AA, JA, IA, UPTR, iFactorA, iFactorB, iFactorC)
       implicit none
       integer, intent (IN)  :: blockID
       integer, intent (IN)  :: NZ
       integer, intent (IN)  :: N 
       real,    intent (IN)  :: dt
       real,    intent (IN)  :: theta  
       real,    intent (OUT) :: AA  (NZ)
       integer, intent (OUT) :: JA  (NZ)
       integer, intent (OUT) :: IA  (N+1)
       integer, intent (OUT) :: UPTR(N)  
       integer, intent(IN)   :: iFactorA
       integer, intent(IN)   :: iFactorB
       integer, OPTIONAL, intent(IN) :: iFactorC 
     end subroutine diff_computeblkAMat
  end interface


end Module diff_interface
