!!****if* source/physics/IncompNS/localAPI/ins_interface
!!
!! NAME
!!  ins_inteface
!!
!! SYNOPSIS
!!  use ins_interface
!!
!! DESCRIPTION
!!  This is an interface module that defines explicit interfaces for private procedures
!!  within the Incompressible Navier Stokes unit.
!!
!!***
Module ins_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"


  ! Routines:
  ! ins_computeDtLocal
  ! ins_ab2rk3.F90
  ! ins_advance.F90
  ! ins_fluxes.F90
  ! ins_rhs.F90
  ! ins_vt.F90


  interface 

    subroutine ins_computeDtLocal(blockID,   & 
                          isize, jsize, ksize,  &
                          dx, dy, dz,           &
                          blkLimits,blkLimitsGC,&
                          facexData,faceyData,  &
                          facezData,            &
                          dtLocal, lminloc)
      implicit none
      integer, intent(IN) :: blockID
      integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
      integer :: isize,jsize,ksize
      real :: dx, dy, dz
      real, intent(INOUT) :: dtLocal
      integer, intent(INOUT) :: lminloc(5)
      real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
    end subroutine ins_computeDtLocal

  end interface



  interface
     subroutine ins_ab2rk3( blockCount, blockList, timeEndAdv, dt)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(INOUT) ::  blockCount
       integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
       real,    INTENT(IN) :: timeEndAdv,dt
       !! -----------------------------------------------------
     end subroutine ins_ab2rk3
  end interface


  interface

      SUBROUTINE ins_predictor(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)
      implicit none 
      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: uni,vni,wni
      REAL :: gama,rhoa,alfa
      END SUBROUTINE ins_predictor

      SUBROUTINE ins_divergence(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)
      implicit none
      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv
      END SUBROUTINE ins_divergence

      SUBROUTINE ins_corrector(uni,vni,wni,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa)
      implicit none
      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz,alfa
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: uni,vni,wni
      END SUBROUTINE ins_corrector

  end interface


  interface

     SUBROUTINE ins_fluxfix(ng,nxc,nyc,nzc,nxi,nyj,nzk,blockCount,&
                            blockList)
     implicit none
     integer, intent(IN) :: ng,nxc,nyc,nzc,nxi,nyj,nzk, &
                            blockCount     
     integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
     end SUBROUTINE ins_fluxfix

     SUBROUTINE ins_fluxfix_p(ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar, &
                              blockCount,blockList)
     implicit none
     integer, intent(in) :: ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar,&
                            blockCount
     integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
     end SUBROUTINE ins_fluxfix_p

  end interface

  interface

      SUBROUTINE ins_rhs3d(uni,vni,wni,tv,ru1,ix1,ix2,jy1,jy2,kz1,kz2,&
                         dx,dy,dz,ru,rv,rw) 
      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, tv
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv, rw
      end SUBROUTINE ins_rhs3d


      SUBROUTINE ins_rhs2d(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv)
      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv
      end SUBROUTINE ins_rhs2d


  end interface


  interface

     SUBROUTINE ins_vt(isgs,ng,nxc,nyc,nzc,RU1,dx,dy,dz,   &
                    coord,bsize,&
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)             

       implicit none
       integer isgs,ng,nxc,nyc,nzc
       real  RU1,dx,dy,dz       
       real coord(3),bsize(3)

       real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
     END SUBROUTINE ins_vt

  end interface


  interface

      subroutine ins_velgradtensor(ng,facevarxx,facevaryy,facevarzz, &
                 dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,       &
                 tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,tpdwdyc,tpdwdzc)
     
      implicit none
      integer ng     
      real dx,dy,dz
      real, pointer, dimension(:,:,:,:) :: facevarxx,facevaryy,facevarzz

      real, dimension(NXB+2*ng,NYB+2*ng,NZB+2*ng) :: tpdudxc,&
             tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
             tpdwdyc,tpdwdzc

      end subroutine ins_velgradtensor

   end interface

   interface

      subroutine ins_setInterpValsGcell(setval)
      logical, intent(IN) :: setval
      end subroutine ins_setInterpValsGcell

   end interface

   interface

      subroutine ins_computeQinout( blockCount, blockList, inou_flg,  Qinout)
        implicit none
        integer, INTENT(IN) :: blockCount
        integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
        logical, INTENT(IN) :: inou_flg
        real,    INTENT(OUT) :: Qinout
      end subroutine

   end interface

   interface

      subroutine ins_rescaleVelout( blockCount, blockList, Qin,  Qout)
        implicit none
        integer, INTENT(IN) :: blockCount
        integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
        real,    INTENT(IN) :: Qin, Qout
      end subroutine

   end interface

   interface

      subroutine ins_convectVelout( blockCount, blockList, convvel)
        implicit none
        integer, INTENT(IN) :: blockCount
        integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
        real,    INTENT(OUT) :: convvel(LOW:HIGH,MDIM)
      end subroutine

   end interface

   interface
      subroutine ins_UstarStats( blockCount, blockList, print_flg, qin_flg)
        implicit none
        !! ---- Argument List ----------------------------------
        integer, INTENT(IN) :: blockCount
        integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
        logical, INTENT(IN) :: print_flg, qin_flg
        !! -----------------------------------------------------
      end subroutine ins_UstarStats

   end interface

   interface

      subroutine ins_pressgradients(time,dt)
        implicit none
        !! ---- Argument List ----------------------------------
        real, intent(in) :: time, dt
        !! -----------------------------------------------------
      end subroutine ins_pressgradients

   end interface

   
   interface 

      subroutine ins_getBulkVelocity(velB,myaxis)
        implicit none
        real, intent(out) :: velB
        integer, intent(in) :: myaxis
      end subroutine ins_getBulkVelocity

   end interface

   interface
      subroutine ins_statsInit()
        implicit none
      end subroutine Ins_statsInit
   end interface



 end Module ins_interface
