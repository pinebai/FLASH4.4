!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_velgradtensor
!!
!! NAME
!!
!!  ins_velgradtensor
!!
!! SYNOPSIS
!!
!!  call ins_velgradtensor(:: ng,
!!                         real, pointer, dimension(:,:,:,:)  :: facevarxx,
!!                         real, pointer, dimension(:,:,:,:)  :: facevaryy,
!!                         real, pointer, dimension(:,:,:,:)  :: facevarzz,
!!                          :: dx,
!!                          :: dy,
!!                          :: dz,
!!                          :: tpdudxc,
!!                          :: tpdudyc,
!!                          :: tpdudzc,
!!                          :: tpdvdxc,
!!                          :: tpdvdyc,
!!                          :: tpdvdzc,
!!                          :: tpdwdxc,
!!                          :: tpdwdyc,
!!                          :: tpdwdzc)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ng : 
!!
!!   facevarxx : 
!!
!!   facevaryy : 
!!
!!   facevarzz : 
!!
!!   dx : 
!!
!!   dy : 
!!
!!   dz : 
!!
!!   tpdudxc : 
!!
!!   tpdudyc : 
!!
!!   tpdudzc : 
!!
!!   tpdvdxc : 
!!
!!   tpdvdyc : 
!!
!!   tpdvdzc : 
!!
!!   tpdwdxc : 
!!
!!   tpdwdyc : 
!!
!!   tpdwdzc : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Subroutine velocgradtensor
! Subroutine to calculate the components of the velocity gradient
! tensor in a block. The values dui/dxj are calculated in cell centers
! and the first line of guardcells, but not in guardcell edges.
!
! Written by Marcos Vanella, December 2006.
! --------------------------------------------------------------------

      subroutine ins_velgradtensor(ng,facevarxx,facevaryy,facevarzz, &
                 dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,       &
                 tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,tpdwdyc,tpdwdzc)
     
      implicit none
#include "constants.h"
#include "Flash.h"

      integer ng 
      real dx,dy,dz
      real, pointer, dimension(:,:,:,:) :: facevarxx,facevaryy,facevarzz

      real, dimension(NXB+2*ng,NYB+2*ng,NZB+2*ng) :: tpdudxc,&
             tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
             tpdwdyc,tpdwdzc


      ! Local variables:
      real tpdudye(NXB+2*ng+1,NYB+2*ng-1,NZB+2*ng),&
           tpdudze(NXB+2*ng+1,NYB+2*ng,NZB+2*ng-1)

      real tpdvdxe(NXB+2*ng-1,NYB+2*ng+1,NZB+2*ng),&
           tpdvdze(NXB+2*ng,NYB+2*ng+1,NZB+2*ng-1)

      real tpdwdxe(NXB+2*ng-1,NYB+2*ng,NZB+2*ng+1),&
           tpdwdye(NXB+2*ng,NYB+2*ng-1,NZB+2*ng+1)

      integer, parameter :: ivel = VELC_FACE_VAR


      ! Initialize
      tpdudxc = 0.
      tpdudyc = 0.
      tpdudzc = 0.
      tpdvdxc = 0.
      tpdvdyc = 0.
      tpdvdzc = 0.
      tpdwdxc = 0.
      tpdwdyc = 0.
      tpdwdzc = 0.


      ! Velocity derivatives:
      ! -------- -----------

      ! du/dx, dv/dy, dw/dz in centers:
      tpdudxc = dx**(-1.)*(facevarxx(ivel,2:NXB+2*ng+1,:,:) - &
                        facevarxx(ivel,1:NXB+2*ng,:,:))

      tpdvdyc = dy**(-1.)*(facevaryy(ivel,:,2:NYB+2*ng+1,:) - &
                        facevaryy(ivel,:,1:NYB+2*ng,:))

      tpdwdzc = dz**(-1.)*(facevarzz(ivel,:,:,2:NZB+2*ng+1) - &
                        facevarzz(ivel,:,:,1:NZB+2*ng))



      ! du/dy: tpdudye(NXB+2*ng+1,NYB+2*ng-1,NZB+2*ng) in edge
      !        tpdudyc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdudye = dy**(-1.)*(facevarxx(ivel,: ,2:NYB+2*ng   ,:) - &
                        facevarxx(ivel,: ,1:NYB+2*ng-1 ,:))

      tpdudyc(:,2:NYB+2*ng-1,:) =                              &
                0.25*(tpdudye(1:NXB+2*ng   ,1:NYB+2*ng-2 ,:) + &
                      tpdudye(2:NXB+2*ng+1 ,1:NYB+2*ng-2 ,:) + &
                      tpdudye(2:NXB+2*ng+1 ,2:NYB+2*ng-1 ,:) + &
                      tpdudye(1:NXB+2*ng   ,2:NYB+2*ng-1 ,:)) 

      ! du/dz: tpdudze(NXB+2*ng+1,NYB+2*ng,NZB+2*ng-1) in edge
      !        tpdudzc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdudze = dz**(-1.)*(facevarxx(ivel,:,:,2:NZB+2*ng) - &
                        facevarxx(ivel,:,:,1:NZB+2*ng-1))

      tpdudzc(:,:,2:NZB+2*ng-1) =                              &
               0.25*(tpdudze(1:NXB+2*ng   ,: ,1:NZB+2*ng-2 ) + &
                     tpdudze(2:NXB+2*ng+1 ,: ,1:NZB+2*ng-2 ) + &
                     tpdudze(2:NXB+2*ng+1 ,: ,2:NZB+2*ng-1 ) + &
                     tpdudze(1:NXB+2*ng   ,: ,2:NZB+2*ng-1 ))            



      ! dv/dxGrid_getBlkCenterCoords: tpdvdxe(NXB+2*ng-1,NYB+2*ng+1,NZB+2*ng) in edge
      !        tpdvdxc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdvdxe = dx**(-1.)*(facevaryy(ivel,2:NXB+2*ng   ,:,:) -      &
                        facevaryy(ivel,1:NXB+2*ng-1 ,:,:))

      tpdvdxc(2:NXB+2*ng-1,:,:) =                             &
                0.25*(tpdvdxe(1:NXB+2*ng-2 ,1:NYB+2*ng   ,:) +&
                      tpdvdxe(2:NXB+2*ng-1 ,1:NYB+2*ng   ,:) +&
                      tpdvdxe(2:NXB+2*ng-1 ,2:NYB+2*ng+1 ,:) +&
                      tpdvdxe(1:NXB+2*ng-2 ,2:NYB+2*ng+1 ,:))  

      ! dv/dz: tpdvdze(NXB+2*ng,NYB+2*ng+1,NZB+2*ng-1) in edge
      !        tpdvdzc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdvdze = dz**(-1.)*(facevaryy(ivel,: ,: ,2:NZB+2*ng   ) -    &
                           facevaryy(ivel,: ,: ,1:NZB+2*ng-1 ))

      tpdvdzc(:,:,2:NZB+2*ng-1) =                              & 
                0.25*(tpdvdze(: ,1:NYB+2*ng   ,1:NZB+2*ng-2 ) +&
                      tpdvdze(: ,1:NYB+2*ng   ,2:NZB+2*ng-1 ) +&
                      tpdvdze(: ,2:NYB+2*ng+1 ,2:NZB+2*ng-1 ) +&
                      tpdvdze(: ,2:NYB+2*ng+1 ,1:NZB+2*ng-2 ))              


      ! dw/dx: tpdwdxe(NXB+2*ng-1,NYB+2*ng,NZB+2*ng+1) in edge
      !        tpdwdxc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdwdxe = dx**(-1.)*(facevarzz(ivel,2:NXB+2*ng   ,: ,:) -     &
                        facevarzz(ivel,1:NXB+2*ng-1 ,: ,:))

      tpdwdxc(2:NXB+2*ng-1,:,:) =                              &
                0.25*(tpdwdxe(1:NXB+2*ng-2 ,: ,1:NZB+2*ng   ) +&
                      tpdwdxe(2:NXB+2*ng-1 ,: ,1:NZB+2*ng   ) +&
                      tpdwdxe(2:NXB+2*ng-1 ,: ,2:NZB+2*ng+1 ) +&
                      tpdwdxe(1:NXB+2*ng-2 ,: ,2:NZB+2*ng+1 ))  

      ! dw/dy: tpdwdye(NXB+2*ng,NYB+2*ng-1,NZB+2*ng+1) in edge
      !        tpdwdyc(NXB+2*ng,NYB+2*ng,NZB+2*ng)     in center
      tpdwdye = dy**(-1.)*(facevarzz(ivel,: ,2:NYB+2*ng   ,:) -     &
                        facevarzz(ivel,: ,1:NYB+2*ng-1 ,:))

      tpdwdyc(:,2:NYB+2*ng-1,:) =                                &
                0.25*(tpdwdye(:  ,1:NYB+2*ng-2  ,1:NZB+2*ng   ) +&
                      tpdwdye(:  ,1:NYB+2*ng-2  ,2:NZB+2*ng+1 ) +&
                      tpdwdye(:  ,2:NYB+2*ng-1  ,2:NZB+2*ng+1 ) +&
                      tpdwdye(:  ,2:NYB+2*ng-1  ,1:NZB+2*ng   )) 


    End subroutine ins_velgradtensor
