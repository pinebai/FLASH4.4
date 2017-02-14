!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/st_calcAccel
!!
!! NAME
!!
!!  st_calcAccel
!!
!! SYNOPSIS
!!
!!  st_calcAccel(integer,INTENT (in)  :: blockid,
!!               integer,dimension(2,MDIM),INTENT(in)  :: blkLimitsGC,
!!               real(:),intent(IN) :: xcoord,
!!               real(:),intent(IN) :: ycoord,
!!               real(:),intent(IN) :: zcoord)
!!
!! DESCRIPTION
!!
!!   Computes components of the zone-averaged forceitational
!!   acceleration.  This version implements a constant forceitational field.
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!   blkLimitsGC : array holding the upper and lower index limits 
!!                 of an entire block (including GC)
!!
!!   xcoord : real array containg the coordinates along IAXIS
!!
!!   ycoord : real array containg the coordinates along JAXIS
!!
!!   zcoord : real array containg the coordinates along KAXIS
!!
!!
!!
!!***

!!REORDER(4): solnData

subroutine st_calcAccel (blockID,blkLimitsGC,xCoord,yCoord,zCoord)

  use Stir_data, ONLY : st_maxmodes, st_mode, st_nmodes, st_aka, st_akb
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"
  
  !     Parameters

  integer,dimension(2,MDIM),INTENT(in)    ::  blkLimitsGC
  integer,INTENT (in) :: blockID

  real,dimension&
       (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),intent(IN)::xCoord
  real,dimension&
       (blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)),intent(IN)::yCoord
  real,dimension&
       (blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),intent(IN)::zCoord
  !     Local
  
  integer       ::  i, j, k, forcedir, m, sizeX
  integer,dimension(MDIM) :: startingPos
  real :: xx, yy, zz
  real :: internalk
  
  ! pre-compute some mode info
  real,dimension(st_maxmodes,&
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS))::cosxi,sinxi
  real,dimension(st_maxmodes,&
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS))::cosxj,sinxj
  real,dimension(st_maxmodes,&
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))::cosxk,sinxk
  real,pointer,dimension(:,:,:,:) :: solnData

  real    :: realtrigterms, imtrigterms
  integer :: ib, ie
  !         precompute some trig.

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  do m = 1, st_nmodes 
     cosxk(m,:) = cos(st_mode(3,m)*zCoord(:))
     sinxk(m,:) = sin(st_mode(3,m)*zCoord(:))
     
     cosxj(m,:) = cos(st_mode(2,m)*yCoord(:))
     sinxj(m,:) = sin(st_mode(2,m)*yCoord(:))
     
     cosxi(m,:) = cos(st_mode(1,m)*xCoord(:))
     sinxi(m,:) = sin(st_mode(1,m)*xCoord(:))
  enddo

  ib=blkLimitsGC(LOW,IAXIS)
  ie=blkLimitsGC(HIGH,IAXIS)
  call Grid_getBlkPtr(blockID,solnData)
  do k = blkLimitsGC (LOW, KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC (LOW, JAXIS), blkLimitsGC(HIGH,JAXIS)
        solnData(ACCX_VAR,ib:ie,j,k)=0.
        solnData(ACCY_VAR,ib:ie,j,k)=0.
        solnData(ACCZ_VAR,ib:ie,j,k)=0.

        do i = blkLimitsGC (LOW, IAXIS), blkLimitsGC(HIGH,IAXIS)
           
           do m = 1, st_nmodes 
              !
              !  these are the real and imaginary parts, respectively, of
              !     e^{ i \vec{k} \cdot \vec{x} }  
              !          = cos(kx*x + ky*y + kz*z) + i sin(kx*x + ky*y + kz*z)
              !
              realtrigterms = ( cosxi(m,i)*cosxj(m,j)*cosxk(m,k)&
                   - sinxi(m,i)*sinxj(m,j)*cosxk(m,k)  & 
                   -sinxi(m,i)*cosxj(m,j)*sinxk(m,k)  & 
                   -cosxi(m,i)*sinxj(m,j)*sinxk(m,k)  )
              imtrigterms = (  cosxi(m,i)*cosxj(m,j)*sinxk(m,k)  & 
                   + cosxi(m,i)*sinxj(m,j)*cosxk(m,k)  & 
                   + sinxi(m,i)*cosxj(m,j)*cosxk(m,k)  & 
                   - sinxi(m,i)*sinxj(m,j)*sinxk(m,k) )
              
                 
                 solnData(ACCX_VAR,i,j,k)=solnData(ACCX_VAR,i,j,k) + 2.*st_aka(IAXIS,m) * realtrigterms & 
                             - 2.*st_akb(IAXIS,m) * imtrigterms
                 solnData(ACCY_VAR,i,j,k)=solnData(ACCY_VAR,i,j,k) + 2.*st_aka(JAXIS,m) * realtrigterms & 
                             - 2.*st_akb(JAXIS,m) * imtrigterms
                 solnData(ACCZ_VAR,i,j,k)=solnData(ACCZ_VAR,i,j,k) + 2.*st_aka(KAXIS,m) * realtrigterms & 
                             - 2.*st_akb(KAXIS,m) * imtrigterms

           enddo  ! end loop over modes

        enddo  ! end loop over i
     enddo  ! end loop over j
  enddo  ! end loop over k
  call Grid_releaseBlkPtr(blockID,solnData)
  return
end subroutine st_calcAccel
