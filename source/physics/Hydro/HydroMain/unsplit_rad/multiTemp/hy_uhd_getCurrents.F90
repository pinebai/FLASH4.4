!!****if* source/physics/Hydro/HydroMain/unsplit_rad/multiTemp/hy_uhd_getCurrents
!!
!! NAME
!!
!!  hy_uhd_getCurrents
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_getCurrents(integer(in) :: blockID,
!!                          integer(in) :: rangeSwitch,
!!                          integer(in) :: blkLimits  (LOW:HIGH,MDIM),
!!                          integer(in) :: datasize(MDIM),
!!                          real(in)    :: del(MDIM),
!!                          real(inout) :: Jp(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS)),  
!!                          real(inout) :: Jm(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS)),  
!!                          integer(in) :: mode_switch,
!!                          real,POINTER,dimension(:,:,:,:) :: scrch_Ptr,
!!                          integer(in),OPTIONAL :: ix,iy,iz) 
!!
!! DESCRIPTION
!! 
!! This routine loops over i,j,k and calculates the
!! current components as J = curlB in various
!! geometries at the respective cell faces. Uses cell
!! centered differencing for both 8-wave and CT
!! Note that J is first order in time, calculated
!! from fields at t = n and not n+1/2, since the
!! latter are not available
!!
!! Scratch Notation for this file: XN01 = Bx,
!! XN02 = By, XN03 = Bz, XN04 = abar, XN05 = zbar
!! -date modified: 10/9/2012 
!! -date modified: April 2013
!!                 - added code for mode_switch==3
!!                 - added code for mode_switch==4
!!
!! ARGUMENTS
!!
!! blockID         - block ID
!! rangeSwitch    - range switch
!! blkLimits       - block limits
!! datasize        - datasize fo allocations
!! mode_switch     - mode switch of calculation method
!! del             - cell's delta 
!! Jp              - current at i+1/2
!! Jm              - current at i-1/2  
!! ix,iy,iz        - i, j, k where the calculation is to be performed (mode 4)
!!
!! PARAMETERS
!!
!!
!!
!!
!!***

subroutine hy_uhd_getCurrents( &
     blockID, rangeSwitch, blkLimits,datasize, del, Jp, Jm, mode_switch,&
     scrch_Ptr,&
     ix,iy,iz)
     
#include "Flash.h"
#include "constants.h"
#include "UHD.h"
#include "Eos.h"

  use Hydro_data,           ONLY : hy_geometry,  &
                                   hy_qele,      &
                                   hy_avogadro
  
  use Grid_interface,       ONLY : Grid_getBlkPtr, &
                                   Grid_releaseBlkPtr, &
                                   Grid_getBlkData, &
                                   Grid_getCellCoords, &
                                   Grid_getBlkIndexLimits
                                     
  implicit none

  ! Arguments:
  integer,intent(in) :: blockID, rangeSwitch
  integer,intent(in) :: blkLimits  (LOW:HIGH,MDIM)
  integer,intent(in) :: datasize(MDIM), mode_switch
  real,   intent(in) :: del(MDIM)

#ifdef FIXEDBLOCKSIZE
  real, intent(inout) :: Jp(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(inout) :: Jm(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
  real, intent(inout) :: Jp(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
  real, intent(inout) :: Jm(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
#endif
  
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr 
  integer, intent(in), OPTIONAL :: ix,iy,iz 
  

#ifdef FLASH_USM_MHD  
  ! Local variables:
  integer :: i,j,k, kx, ky, kz
  integer :: imin,imax,jmin,jmax,kmin,kmax
  integer :: blklmts(LOW:HIGH,MDIM)
  integer :: blklmtsgc(LOW:HIGH,MDIM)
 
  real, pointer, dimension(:,:,:,:) :: U
  
  real    :: dx, dy, dz
  real, allocatable :: xc(:)
  integer :: isize, jsize, ksize, iskip
  real :: dyBxp, dzBxp, dxByp, dzByp, dxBzp, dyBzp
  real :: dyBxm, dzBxm, dxBym, dzBym, dxBzm, dyBzm
  real :: abar, zbar
  real :: pelep, pionp, pradp
  real :: pelem, pionm, pradm
  
  ! Subroutine body:
  isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1


  !! Set ranges for update
  imin  = blkLimits(LOW, IAXIS)
  imax  = blkLimits(HIGH,IAXIS)
  jmin  = 1
  jmax  = 1
  kmin  = 1
  kmax  = 1

  dx = del(DIR_X)
  dy = 1.
  dz = 1.
  if (NDIM >= 2) then
     jmin  = blkLimits(LOW, JAXIS)
     jmax  = blkLimits(HIGH,JAXIS)
     dy = del(DIR_Y)
     if (NDIM == 3) then
        kmin  = blkLimits(LOW, KAXIS)
        kmax  = blkLimits(HIGH,KAXIS)
        dz = del(DIR_Z)
     endif
  endif

  ! initialize B derivatives to zero
  dxByp = 0.0
  dxBzp = 0.0 
  dyBxp = 0.0 
  dyBzp = 0.0
  dzBxp = 0.0
  dzByp = 0.0 

  dxBym = 0.0
  dxBzm = 0.0 
  dyBxm = 0.0 
  dyBzm = 0.0
  dzBxm = 0.0
  dzBym = 0.0 
  
  ! get limits 
  call Grid_getBlkIndexLimits(blockId,blklmts,blklmtsgc)
  
  ! allocate coords
  allocate(xc(blklmtsgc(HIGH, IAXIS)))
  
  ! get coords
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xc, blklmtsgc(HIGH, IAXIS))
  

  ! get block pointers
  call Grid_getBlkPtr(blockID,U,CENTER)
  
  ! define dimension dependent switches
  kx=0
  ky=0
  kz=0

  if (NDIM > 1) then
     ky=1
     if (NDIM > 2) then
        kz=1
     endif
  endif

  iskip = 1
  if (NDIM == 1 .and. rangeSwitch .eq. UPDATE_BOUND) iskip = imax-imin


  if (mode_switch < 4) then
     do k=kmin-kx*kz,kmax+kx*kz
        do j=jmin-kx*ky,jmax+kx*ky
           if (NDIM >= 2) then
              iskip = 1
              if (rangeSwitch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
                 iskip = imax-imin
                 if (NDIM == 3) then
                    iskip = 1
                    if (k > kmin .and. k < kmax) then
                       iskip = imax-imin
                    endif
                 endif
              endif
           endif
           do i=imin-kx,imax+kx,iskip
      
           ! define abar & zbar from the scratch array we saved in unsplitUpdate
           abar= scrch_Ptr(HY_XN04_SCRATCHCTR_VAR,i,j,k)
           zbar=scrch_Ptr(HY_XN05_SCRATCHCTR_VAR,i,j,k) 
        
           if (hy_geometry == CARTESIAN) then

  !* ********************************************* *
  !*   calculate currents as 
  !*   Jp(1,i,j,k) = dyBzp - dzByp @ i+1/2
  !*   Jp(2,i,j,k) = dzBxp - dxBzp @ j+1/2
  !*   Jp(3,i,j,k) = dxByp - dyBxp @ k+1/2    
  !*
  !*   Jm(1,i,j,k) = dyBzm - dzBym @ i-1/2
  !*   Jm(2,i,j,k) = dzBxm - dxBzm @ j-1/2
  !*   Jm(3,i,j,k) = dxBym - dyBxm @ k-1/2    
  !* ********************************************* *        

              if (NDIM == 3) dxByp = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k+1) &
                                   -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k+1) &
                                   +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k  ) &
                                   -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k  ))/dx !@ k+1/2
              if (NDIM > 1)  dxBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j+1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j+1,k) &
                                   +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j  ,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j  ,k))/dx !@ j+1/2 

              if (NDIM == 3) dxBym = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k-1) &
                               -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k-1) &
                               +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k  ) &
                               -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k  ))/dx !@ k-1/2
              if (NDIM > 1)  dxBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j-1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j-1,k) &
                                   +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j  ,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j  ,k))/dx !@ j-1/2 

              if (NDIM > 1) then

                 if (NDIM == 3) dyBxp = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k+1) &
                                      -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k+1) &
                                      +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k  ) &
                                      -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k  ))/dy !@ k+1/2
                 dyBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j+1,k) &
                       -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j-1,k) &  
                       +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j+1,k) &
                       -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j-1,k))/dy !@ i+1/2

                 if (NDIM == 3) dyBxm = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k-1) &
                                      -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k-1) &
                                      +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k  ) &
                                      -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k  ))/dy !@ k-1/2
                 dyBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j+1,k) &
                       -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j-1,k) &  
                       +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j+1,k) &
                       -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j-1,k))/dy !@ i-1/2
            
              endif
              if (NDIM == 3) then
            
                 dzBxp = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k+1) &
                       -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j+1,k-1) &  
                       +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j  ,k+1) &
                       -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j  ,k-1))/dz !@ j+1/2
            
                 dzByp = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k+1) &
                       -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i+1,j,k-1) &  
                       +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i  ,j,k+1) &
                       -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i  ,j,k-1))/dz !@ i+1/2

                 dzBxm = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k+1) &
                       -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j-1,k-1) &  
                       +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j  ,k+1) &
                       -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,i,j  ,k-1))/dz !@ j-1/2
            
                 dzBym = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k+1) &
                       -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i-1,j,k-1) &  
                       +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i  ,j,k+1) &
                       -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,i  ,j,k-1))/dz !@ i-1/2
              endif

               Jp(1,i,j,k) = dyBzp - dzByp !@ i+1/2
               Jp(2,i,j,k) = dzBxp - dxBzp !@ j+1/2
               Jp(3,i,j,k) = dxByp - dyBxp !@ k+1/2    
  
               Jm(1,i,j,k) = dyBzm - dzBym !@ i-1/2
               Jm(2,i,j,k) = dzBxm - dxBzm !@ j-1/2
               Jm(3,i,j,k) = dxBym - dyBxm !@ k-1/2        
           endif   

           if (hy_geometry == CYLINDRICAL) then
  !* ********************************************* *
  !*   In this geometry (2D) we require only J1 and J2 for
  !*   the R and Z components, which are calculated as
  !*   (** note: x=R, y=Z, z=phi **)
  !*
  !*   Jp(1,i,j,k) = -dyBzp = - dyBzp @ i+1/2
  !*   Jp(2,i,j,k) =  dxBzp = (1/xc) dx(xc Bz)p @ j+1/2
  !*
  !*   Jm(1,i,j,k) = -dyBzm = - dyBzm @ i-1/2
  !*   Jm(2,i,j,k) =  dxBzm = (1/xc) dx(xc Bz)m @ j-1/2
  !* ********************************************* *        
              if (NDIM > 1)  dyBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j+1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j-1,k) &
                                   +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j+1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j-1,k))/dy !@ i+1/2 
                              
              if (NDIM > 1)  dxBzp = 0.25*(1./xc(i))*(xc(i+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j+1,k) &
                                   -                  xc(i-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j+1,k) &
                                   +                  xc(i+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j  ,k) &
                                   -                  xc(i-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j  ,k))/dx !@ j+1/2 

              if (NDIM > 1)  dyBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j+1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j-1,k) &
                                   +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j+1,k) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i  ,j-1,k))/dy !@ i-1/2 
                               
              if (NDIM > 1)  dxBzm = 0.25*(1./xc(i))*(xc(i+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j-1,k) &
                                   -                  xc(i-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j-1,k) &
                                   +                  xc(i+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i+1,j  ,k) &
                                   -                  xc(i-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,i-1,j  ,k))/dx !@ j-1/2 
            
        
              Jp(1,i,j,k) = -dyBzp !@ i+1/2
              Jp(2,i,j,k) =  dzBxp !@ j+1/2
  
              Jm(1,i,j,k) = -dyBzm !@ i-1/2
              Jm(2,i,j,k) =  dzBxm !@ j-1/2
          
           endif   
      
  !* ********************************************* *
  !*   Having calculated J, we now multiply it with
  !*   Eele and divide by (ne qe) to obtain the
  !*   correction term for the fluxes of Eele in
  !*   unsplitUpdateMultiTemp. Comment this part
  !*   or introduce a switch if you want to use
  !*   this routine for normal J calculation.
  !* ********************************************* *        

  ! multiply and divide with average values
  ! Note: U(EELE_VAR) is still the t=n value
  ! since the update occurs at unsplitUpdateMultiTemp
  ! BUT density has already been updated so we get it
  ! from scratch, VAR1 
           if (mode_switch == 1) then ! UPDATE J TO INCLUDE CORRECTION OF E_ELE
              Jp(1,i,j,k) = Jp(1,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i+1,j,k))
!             Jp(1,i,j,k) = Jp(1,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i+1,j,k)))
              Jp(1,i,j,k) = Jp(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
        
              Jm(1,i,j,k) = Jm(1,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i-1,j,k))
!             Jm(1,i,j,k) = Jm(1,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i-1,j,k)))
              Jm(1,i,j,k) = Jm(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                
              if (NDIM >1) then

                 Jp(2,i,j,k) = Jp(2,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j+1,k))
!                Jp(2,i,j,k) = Jp(2,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j+1,k)))
                 Jp(2,i,j,k) = Jp(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
        
                 Jm(2,i,j,k) = Jm(2,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j-1,k))
!                Jm(2,i,j,k) = Jm(2,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j-1,k)))
                 Jm(2,i,j,k) = Jm(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
        
                 if (NDIM ==3) then

                    Jp(3,i,j,k) = Jp(3,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j,k+1))
!                   Jp(3,i,j,k) = Jp(3,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k+1)))
                    Jp(3,i,j,k) = Jp(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
        
                    Jm(3,i,j,k) = Jm(3,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j,k-1))
!                   Jm(3,i,j,k) = Jm(3,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k-1)))
                    Jm(3,i,j,k) = Jm(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
          
                 endif
              endif 
           else if (mode_switch==2) then! UPDATE J TO INCLUDE CORRECTION OF E_TOT
        
            ! recover pelem pelep
              call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                       U(:,i,j,k), pelem, pionm, pradm)
              call hy_uhd_getPressure( U(EELE_VAR,i+1,j,k), U(EION_VAR,i+1,j,k), U(ERAD_VAR,i+1,j,k), &
                                       U(:,i+1,j,k), pelep, pionp, pradp)

        
              Jp(1,i,j,k) = Jp(1,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i+1,j,k) + pelep + pelem)
!             Jp(1,i,j,k) = Jp(1,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i+1,j,k)))
              Jp(1,i,j,k) = Jp(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
              call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                       U(:,i,j,k), pelep, pionp, pradp)
              call hy_uhd_getPressure( U(EELE_VAR,i-1,j,k), U(EION_VAR,i-1,j,k), U(ERAD_VAR,i-1,j,k), &
                                        U(:,i-1,j,k), pelem, pionm, pradm)
        
              Jm(1,i,j,k) = Jm(1,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i-1,j,k) + pelep + pelem)
!             Jm(1,i,j,k) = Jm(1,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i-1,j,k)))
              Jm(1,i,j,k) = Jm(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                
              if (NDIM >1) then
    
            ! recover pelem pelep
                 call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                           U(:,i,j,k), pelem, pionm, pradm)
                 call hy_uhd_getPressure( U(EELE_VAR,i,j+1,k), U(EION_VAR,i,j+1,k), U(ERAD_VAR,i,j+1,k), &
                                           U(:,i,j+1,k), pelep, pionp, pradp)

                 Jp(2,i,j,k) = Jp(2,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j+1,k) + pelep + pelem)
!                Jp(2,i,j,k) = Jp(2,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j+1,k)))
                 Jp(2,i,j,k) = Jp(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
                 call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                     U(:,i,j,k), pelep, pionp, pradp)
                 call hy_uhd_getPressure( U(EELE_VAR,i,j-1,k), U(EION_VAR,i,j-1,k), U(ERAD_VAR,i,j-1,k), &
                                      U(:,i,j-1,k), pelem, pionm, pradm)
        
                 Jm(2,i,j,k) = Jm(2,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j-1,k) + pelep + pelem)
!                Jm(2,i,j,k) = Jm(2,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j-1,k)))
                 Jm(2,i,j,k) = Jm(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
        
                    if (NDIM ==3) then

              ! recover pelem pelep
                       call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                                U(:,i,j,k), pelem, pionm, pradm)
                       call hy_uhd_getPressure( U(EELE_VAR,i,j,k+1), U(EION_VAR,i,j,k+1), U(ERAD_VAR,i,j,k+1), &
                                                U(:,i,j,k+1), pelep, pionp, pradp)

                       Jp(3,i,j,k) = Jp(3,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j,k+1) + pelep + pelem)
!                      Jp(3,i,j,k) = Jp(3,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k+1)))
                       Jp(3,i,j,k) = Jp(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
                       call hy_uhd_getPressure( U(EELE_VAR,i,j,k), U(EION_VAR,i,j,k), U(ERAD_VAR,i,j,k), &
                                                U(:,i,j,k), pelep, pionp, pradp)
                       call hy_uhd_getPressure( U(EELE_VAR,i,j,k-1), U(EION_VAR,i,j,k-1), U(ERAD_VAR,i,j,k-1), &
                                                U(:,i,j,k-1), pelem, pionm, pradm)
        
                       Jm(3,i,j,k) = Jm(3,i,j,k)*0.5*(U(EELE_VAR,i,j,k) +U(EELE_VAR,i,j,k-1) + pelep + pelem)
!                      Jm(3,i,j,k) = Jm(3,i,j,k)/(0.5*(scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)+scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k-1)))
                       Jm(3,i,j,k) = Jm(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                    endif
                 endif
              else if (mode_switch==3) then! UPDATE J TO INCLUDE CORRECTION OF V 
        
        
                 Jp(1,i,j,k) = Jp(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                 Jm(1,i,j,k) = Jm(1,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
               
                 if (NDIM >1) then
   
                    Jp(2,i,j,k) = Jp(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                    Jm(2,i,j,k) = Jm(2,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
       
                    if (NDIM ==3) then

                       Jp(3,i,j,k) = Jp(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                       Jm(3,i,j,k) = Jm(3,i,j,k)*abar/(zbar*hy_avogadro*hy_qele) 
                    endif
                 endif
              endif
           end do
        end do
     end do
  endif
  if (mode_switch == 4) then !! case 4 is for the update of total energy, given ix,iy,iz when called within a loop!  

     abar= scrch_Ptr(HY_XN04_SCRATCHCTR_VAR,ix,iy,iz)
     zbar=scrch_Ptr(HY_XN05_SCRATCHCTR_VAR,ix,iy,iz) 
        
     if (hy_geometry == CARTESIAN) then

  !* ********************************************* *
  !*   calculate currents as 
  !*   Jp(1,i,j,k) = dyBzp - dzByp @ i+1/2
  !*   Jp(2,i,j,k) = dzBxp - dxBzp @ j+1/2
  !*   Jp(3,i,j,k) = dxByp - dyBxp @ k+1/2    
  !*
  !*   Jm(1,i,j,k) = dyBzm - dzBym @ i-1/2
  !*   Jm(2,i,j,k) = dzBxm - dxBzm @ j-1/2
  !*   Jm(3,i,j,k) = dxBym - dyBxm @ k-1/2    
  !* ********************************************* *        

        if (NDIM == 3) dxByp = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz+1) &
                             -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz+1) &
                             +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz  ) &
                             -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz  ))/dx !@ k+1/2
        if (NDIM > 1)  dxBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy+1,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy+1,iz) &
                             +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy  ,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy  ,iz))/dx !@ j+1/2 

        if (NDIM == 3) dxBym = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz-1) &
                             -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz-1) &
                             +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz  ) &
                             -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz  ))/dx !@ k-1/2
        if (NDIM > 1)  dxBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy-1,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy-1,iz) &
                             +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy  ,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy  ,iz))/dx !@ j-1/2 

        if (NDIM > 1) then

           if (NDIM == 3) dyBxp = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz+1) &
                                -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz+1) &
                                +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz  ) &
                                -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz  ))/dy !@ k+1/2
           dyBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy+1,iz) &
                 -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy-1,iz) &  
                 +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy+1,iz) &
                 -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy-1,iz))/dy !@ i+1/2

           if (NDIM == 3) dyBxm = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz-1) &
                                -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz-1) &
                                +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz  ) &
                                -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz  ))/dy !@ k-1/2
           dyBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy+1,iz) &
                 -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy-1,iz) &  
                 +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy+1,iz) &
                 -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy-1,iz))/dy !@ i-1/2
            
        endif
        if (NDIM == 3) then
            
           dzBxp = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz+1) &
                 -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy+1,iz-1) &  
                 +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy  ,iz+1) &
                 -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy  ,iz-1))/dz !@ j+1/2
            
           dzByp = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz+1) &
                 -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix+1,iy,iz-1) &  
                 +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix  ,iy,iz+1) &
                 -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix  ,iy,iz-1))/dz !@ i+1/2

           dzBxm = 0.25*(scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz+1) &
                 -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy-1,iz-1) &  
                 +       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy  ,iz+1) &
                 -       scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,ix,iy  ,iz-1))/dz !@ j-1/2
            
           dzBym = 0.25*(scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz+1) &
                 -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix-1,iy,iz-1) &  
                 +       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix  ,iy,iz+1) &
                 -       scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,ix  ,iy,iz-1))/dz !@ i-1/2
        endif

        Jp(1,ix,iy,iz) = dyBzp - dzByp !@ i+1/2
        Jp(2,ix,iy,iz) = dzBxp - dxBzp !@ j+1/2
        Jp(3,ix,iy,iz) = dxByp - dyBxp !@ k+1/2    
  
        Jm(1,ix,iy,iz) = dyBzm - dzBym !@ i-1/2
        Jm(2,ix,iy,iz) = dzBxm - dxBzm !@ j-1/2
        Jm(3,ix,iy,iz) = dxBym - dyBxm !@ k-1/2        
     endif   

     if (hy_geometry == CYLINDRICAL) then
  !* ********************************************* *
  !*   In this geometry (2D) we require only J1 and J2 for
  !*   the R and Z components, which are calculated as
  !*   (** note: x=R, y=Z, z=phi **)
  !*
  !*   Jp(1,i,j,k) = -dyBzp = - dyBzp @ i+1/2
  !*   Jp(2,i,j,k) =  dxBzp = (1/xc) dx(xc Bz)p @ j+1/2
  !*
  !*   Jm(1,i,j,k) = -dyBzm = - dyBzm @ i-1/2
  !*   Jm(2,i,j,k) =  dxBzm = (1/xc) dx(xc Bz)m @ j-1/2
  !* ********************************************* *        
        if (NDIM > 1)  dyBzp = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy+1,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy-1,iz) &
                             +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy+1,iz) &
                             -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy-1,iz))/dy !@ i+1/2 
                              
        if (NDIM > 1)  dxBzp = 0.25*(1./xc(ix))*(xc(ix+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy+1,iz) &
                                   -                  xc(ix-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy+1,iz) &
                                   +                  xc(ix+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy  ,iz) &
                                   -                  xc(ix-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy  ,iz))/dx !@ j+1/2 

        if (NDIM > 1)  dyBzm = 0.25*(scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy+1,iz) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy-1,iz) &
                                   +       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy+1,iz) &
                                   -       scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix  ,iy-1,iz))/dy !@ i-1/2 
                               
        if (NDIM > 1)  dxBzm = 0.25*(1./xc(ix))*(xc(ix+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy-1,iz) &
                                   -                  xc(ix-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy-1,iz) &
                                   +                  xc(ix+1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix+1,iy  ,iz) &
                                   -                  xc(ix-1)*scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,ix-1,iy  ,iz))/dx !@ j-1/2 
            
        
        Jp(1,ix,iy,iz) = -dyBzp !@ i+1/2
        Jp(2,ix,iy,iz) =  dzBxp !@ j+1/2
  
        Jm(1,ix,iy,iz) = -dyBzm !@ i-1/2
        Jm(2,ix,iy,iz) =  dzBxm !@ j-1/2
          
     endif   
      
  !* ********************************************* *
  !*   Having calculated J, we now multiply it with
  !*   Eele and divide by (ne qe) to obtain the
  !*   correction term for the fluxes of Eele in
  !*   unsplitUpdateMultiTemp. Comment this part
  !*   or introduce a switch if you want to use
  !*   this routine for normal J calculation.
  !* ********************************************* *        

  ! multiply and divide with average values
  ! Note: U(EELE_VAR) is still the t=n value
  ! since the update occurs at unsplitUpdateMultiTemp
  ! BUT density has already been updated so we get it
  ! from scratch, VAR1 
        
            ! recover pelem pelep
     call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                              U(:,ix,iy,iz), pelem, pionm, pradm)
     call hy_uhd_getPressure( U(EELE_VAR,ix+1,iy,iz), U(EION_VAR,ix+1,iy,iz), U(ERAD_VAR,ix+1,iy,iz), &
                              U(:,ix+1,iy,iz), pelep, pionp, pradp)

        
     Jp(1,ix,iy,iz) = Jp(1,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix+1,iy,iz) + pelep + pelem)
     Jp(1,ix,iy,iz) = Jp(1,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
     call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                              U(:,ix,iy,iz), pelep, pionp, pradp)
     call hy_uhd_getPressure( U(EELE_VAR,ix-1,iy,iz), U(EION_VAR,ix-1,iy,iz), U(ERAD_VAR,ix-1,iy,iz), &
                              U(:,ix-1,iy,iz), pelem, pionm, pradm)
        
     Jm(1,ix,iy,iz) = Jm(1,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix-1,iy,iz) + pelep + pelem)
     Jm(1,ix,iy,iz) = Jm(1,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 
                
     if (NDIM >1) then
    
            ! recover pelem pelep
        call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                                 U(:,ix,iy,iz), pelem, pionm, pradm)
        call hy_uhd_getPressure( U(EELE_VAR,ix,iy+1,iz), U(EION_VAR,ix,iy+1,iz), U(ERAD_VAR,ix,iy+1,iz), &
                                 U(:,ix,iy+1,iz), pelep, pionp, pradp)

        Jp(2,ix,iy,iz) = Jp(2,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix,iy+1,iz) + pelep + pelem)
        Jp(2,ix,iy,iz) = Jp(2,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
        call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                                 U(:,ix,iy,iz), pelep, pionp, pradp)
        call hy_uhd_getPressure( U(EELE_VAR,ix,iy-1,iz), U(EION_VAR,ix,iy-1,iz), U(ERAD_VAR,ix,iy-1,iz), &
                                 U(:,ix,iy-1,iz), pelem, pionm, pradm)
        
        Jm(2,ix,iy,iz) = Jm(2,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix,iy-1,iz) + pelep + pelem)
        Jm(2,ix,iy,iz) = Jm(2,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 
        
        if (NDIM ==3) then

              ! recover pelem pelep
           call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                                    U(:,ix,iy,iz), pelem, pionm, pradm)
           call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz+1), U(EION_VAR,ix,iy,iz+1), U(ERAD_VAR,ix,iy,iz+1), &
                                    U(:,ix,iy,iz+1), pelep, pionp, pradp)

           Jp(3,ix,iy,iz) = Jp(3,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix,iy,iz+1) + pelep + pelem)
           Jp(3,ix,iy,iz) = Jp(3,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 

            ! recover pelem pelep
           call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz), U(EION_VAR,ix,iy,iz), U(ERAD_VAR,ix,iy,iz), &
                                    U(:,ix,iy,iz), pelep, pionp, pradp)
           call hy_uhd_getPressure( U(EELE_VAR,ix,iy,iz-1), U(EION_VAR,ix,iy,iz-1), U(ERAD_VAR,ix,iy,iz-1), &
                                    U(:,ix,iy,iz-1), pelem, pionm, pradm)
        
           Jm(3,ix,iy,iz) = Jm(3,ix,iy,iz)*0.5*(U(EELE_VAR,ix,iy,iz) +U(EELE_VAR,ix,iy,iz-1) + pelep + pelem)
           Jm(3,ix,iy,iz) = Jm(3,ix,iy,iz)*abar/(zbar*hy_avogadro*hy_qele) 
        endif
     endif
  endif
 
  ! Release block pointers
  call Grid_releaseBlkPtr(blockID,U,CENTER)
     
  deallocate(xc)

#endif 
end subroutine hy_uhd_getCurrents
