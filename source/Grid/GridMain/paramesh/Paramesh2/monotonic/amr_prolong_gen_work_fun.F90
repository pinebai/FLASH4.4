

subroutine amr_prolong_gen_work_fun & 
     &       ( ia, ib, ja, jb, ka, kb, isg, & 
     &        ioff, joff, koff, mype)
  

  use workspace
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_data,ONLY: gr_convertToConsvdForMeshCalls, gr_dirGeom, gr_smallx, gr_intpol
  implicit none
#include "constants.h"
#include "Flash.h"
  
  integer,parameter :: niver = 1

  integer,parameter :: lw = 4*(NXB+2*NGUARD+4) + K2D*(4*(NYB+2*NGUARD +4)+&
       5*(NXB+2*NGUARD+4)*NUNK_VARS) + K3D*(4*(NZB+2*NGUARD+4)+&
       5*(NXB+2*NGUARD+4)*(NYB+2*NGUARD+4)*NUNK_VARS)
  integer,parameter :: liw=  NXB+NYB+NZB+12+2*NGUARD


  integer,parameter :: ref_ratio = 2
  integer,parameter :: intpol_guard = 2
  integer,parameter :: twice_iguard = 2*intpol_guard

  integer,parameter :: lxiend = (NXB+2*NGUARD)/2 + twice_iguard
  integer,parameter :: linxu  = lxiend*niver
  integer,parameter :: lyiend = max(1,K2D*((NYB+2*NGUARD)/2 + twice_iguard))
  integer,parameter :: lziend = max(1,K3D*((NZB+2*NGUARD)/2 + twice_iguard))

  integer,parameter :: lxoend = NXB+2*NGUARD
  integer,parameter :: lonxu  = lxoend*niver
  integer,parameter :: lyoend = NYB+2*NGUARD
  integer,parameter :: lzoend = NZB+2*NGUARD

  integer, intent(IN) :: isg, mype, ia, ib, ja, jb, ka, kb,ioff,joff,koff
  
  integer :: i,j,k,ivar,i1,j1,k1,offia,offib,offic,offoa,offob,offoc
 
  
  real, dimension(1:2*NGUARD+NXB)     :: xo, xi
  real, dimension(1:2*NGUARD*K2D+NYB) :: yo, yi
  real, dimension(1:2*NGUARD*K3D+NZB) :: zo, zi

  real, dimension(linxu,lyiend,lziend) :: pu
  real, dimension(lonxu,lyoend,lzoend) :: qu

  real, dimension(lw)  :: wmap    ! the two work arrays needed by the inter-
  integer, dimension(liw) :: iwmap   ! polation routine
  ! note variable above used to real, fixed in paramesh2 and paramesh3 as well
  
  real :: pdx, pdy, pdz, cdx, cdy, cdz ! half grid spacings for parent and child
  
  integer       :: imapcx, ierr
  integer       :: geom(3)
  integer       :: xoend,yoend,zoend,xiend,yiend,ziend
  integer       :: inxu, onxu


#ifdef DEBUG_GRID

  if(NGUARD<twice_iguard) then
     call abort_flash("[AMR_PROLONG_GEN_WORK_FUN] ERROR: not enough guardcells to support asked for interpolation")
  end if
#endif     
  geom=gr_dirGeom

!!  First initialize all values that may not be used in non 3D case
  pdy   = 0.e0
  cdy   = 0.e0
  pdz   = 0.e0
  cdz   = 0.e0
  zoend = 1
  yoend = 1
  ziend = 1
  yiend = 1
  offob = 0
  offoc = 0
  offib = 0
  offic = 0
  
  offia = ia/2 + NGUARD - twice_iguard
  offoa = ia - 1
  
  xoend = ib-ia+1
  xiend = xoend/ref_ratio + twice_iguard  

  call Grid_getCellCoords(IAXIS,isg,CENTER,.true.,xi,GRID_IHI_GC)
!!  call dBaseGetCoords(izn,  iXcoord, isg, xi)
  pdx = xi(2) - xi(1)
  cdx = 0.5e0*pdx

  do i = 1,xoend
     xo(i)=xi(offoa+i)
  end do
  
  xi(1) = xo(1)+ cdx - twice_iguard*pdx 

  do i = 2,xiend
     xi(i) = xi(i-1) + 2.e0*pdx
  end do
#if N_DIM > 1

  offib = ja/2 + NGUARD - twice_iguard
  offob = ja - 1
  yoend = jb-ja+1
  yiend = yoend/ref_ratio + twice_iguard  

  call Grid_getCellCoords(JAXIS,isg,CENTER,.true.,yi,GRID_JHI_GC)
!!  call dBaseGetCoords(izn,  iYcoord, isg, yi)
  pdy = yi(2) - yi(1)
  cdy = 0.5e0*pdy
  
  do i = 1,yoend
     yo(i)=yi(offob+i)
  end do
  yi(1) = yo(1)+cdy -twice_iguard*pdy
  
  do i = 2,yiend
     yi(i) = yi(i-1) + 2.e0*pdy
  end do
#endif
#if N_DIM == 3

  offic = ka/2 + NGUARD - twice_iguard
  offoc = ka-1
  
  zoend = kb-ka+1
  ziend = zoend/ref_ratio + twice_iguard  

  call Grid_getCellCoords(KAXIS,isg,CENTER,.true.,zi,GRID_KHI_GC)
!!  call dBaseGetCoords(izn,  iZcoord, isg, zi)
  pdz = zi(2) - zi(1)
  cdz = 0.5e0*pdz
  
  do i = 1,zoend
     zo(i)=zi(offoc+i)
  end do
  zi(1) = zo(1)+cdz - twice_iguard*pdz
  
  do i = 2,ziend
     zi(i) = zi(i-1) + 2.e0*pdz
  end do
#endif

  pu = 0.e0

  do k = 1,ziend
     k1 = (offic+koff)*K3D+k
     do j = 1,yiend
        j1 = (offib+joff)*K2D+j
        do i = 1,xiend
           i1 = offia+ioff+i
           pu(i,j,k) = recv1(i1,j1,k1)
        end do
     end do
  end do

  inxu = niver*xiend
  onxu = niver*xoend  
     
  imapcx = 1

  if(ndim == 1) then
     call umap1 (xiend,xi,pdx,inxu,pu,&
          &      xoend,xo,cdx,onxu,qu,&
          &      niver,gr_intpol,imapcx,&  
          &      geom(1),ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,gr_convertToConsvdForMeshCalls,gr_smallx)
  elseif(ndim==2) then
     call umap2 (lxiend,lyiend,xiend,yiend,linxu,inxu,xi,pdx,yi,pdy,pu,&
          &      lxoend,lyoend,xoend,yoend,lonxu,onxu,xo,cdx,yo,cdy,qu,&
          &      niver,gr_intpol,imapcx,&  
          &      geom(1),geom(2),ref_ratio,ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,gr_convertToConsvdForMeshCalls,gr_smallx)
  else
     call umap3 (lxiend,lyiend,lziend,xiend,yiend,ziend,&
          &      linxu,inxu,xi,pdx,yi,pdy,zi,pdz,pu,&
          &      lxoend,lyoend,lzoend,xoend,yoend,zoend,&
          &      lonxu,onxu,xo,cdx,yo,cdy,zo,cdz,qu,&
          &      niver,gr_intpol,imapcx,&
          &      geom(1),geom(2),geom(3),ref_ratio,ref_ratio,ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,gr_convertToConsvdForMeshCalls,gr_smallx)
  end if
  
  do k = 1,zoend
     k1 = offoc*K3D+k
     do j = 1,yoend
        j1 = offob*K2D+j
        do ivar = 1,niver
           do i = 1,xoend
              i1 = offoa+i
              work(i1,j1,k1,isg,1) = qu(i,j,k)
           end do
        end do
     end do
  end do

  return
end subroutine amr_prolong_gen_work_fun
