!!****if* source/Simulation/SimulationMain/SBlast/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(  integer(in) :: blockID,
!!                         
!!
!! DESCRIPTION   
!!     Initialize fluid properties (density, pressure, velocity, etc.) 
!!          in a single block 
!!          for the Two gamma problem. There are two fluids in the 
!!          problem that propagate together.
!!
!! ARGUMENTS
!!      blockID    -      the current block number to be filled
!!       
!!
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "Eos.h"


subroutine calc_atmos(atmos, rho_atmos, p_atmos, gamma, sh, r, rho, p)
   use Driver_interface, ONLY : Driver_abortFlash

! Compute density and pressure for a specified type of atmosphere

   implicit none
   integer, intent(IN) :: atmos                           ! tyep of atmosphere
   real, intent(IN) :: rho_atmos, p_atmos, gamma, sh, r   ! input density, pressure, gamma, scale height, and height
   real, intent(OUT) :: rho, p                            ! output density and pressure

   real :: eint, beta

   if (atmos==CONST_RHO_ATMOS) then
      rho = rho_atmos
      p   = p_atmos
   else if (atmos==CONST_P_ATMOS) then
      rho = rho_atmos*exp(-r/sh)
      p   = p_atmos
   else
      call Driver_abortFlash("Invalid value for atmos in calc_atmos")
   endif

end subroutine calc_atmos



subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_putPointData, &
    Grid_getGeometry, Grid_getPointData

  implicit none

  integer, intent(IN) :: blockId

  integer :: i, j, k, n     !loop counters
  
  integer :: istat

  real, allocatable, DIMENSION(:) :: xCenter, yCenter, zCenter
  real, allocatable, DIMENSION(:) :: xLeft, yLeft, zLeft
  real, allocatable, DIMENSION(:) :: xRight, yRight, zRight

  real, pointer, dimension(:,:,:,:):: solnData

  integer, dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  logical :: gcell=.true.
  integer,dimension(MDIM):: axis

  integer :: geo                           

  ! Distance from origin for the center of a cell and each of the cell corners
  real :: r, r_bne, r_bnw, r_bsw, r_bse, r_tne, r_tnw, r_tsw, r_tse

  real :: rho, p, eint, temp, gamc, game          ! Thermodynamic variables density, pressure, temperature, gamc, game
  real :: mfrac1, mfrac2, mfrac3                  ! Mass fractions

  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: massFrac
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask  

  integer :: ii, jj, kk                            ! Subcycle loop counters
  real dx, dy, dz, dr, xx, yy, zz, rr, r2, r1      ! Volume elements, volume centers, and edges for subscycle computations
  real theta1, theta2

  real :: dVol, pmassIn, pmassOut                  ! Volume and mass fractions
  real :: Vol                                      ! Volume

  real :: rho_temp, p_temp, p_sum                  ! Scalar variables used in subcycle computation

  do kk = EOS_VARS+1,EOS_NUM
     mask(kk) = .FALSE.
  enddo

  call Grid_getGeometry(geo)

  ! Initialize unk variable storage to zero:
  !
  !call Grid_getBlkPtr(blockID,solnData)
  !solndata = 0.0
  !call Grid_releaseBlkPtr(solnData)
    
  ! set coordinates
  !
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCenter(sizeX),stat=istat)
  allocate(yCenter(sizeY),stat=istat)
  allocate(zCenter(sizeZ),stat=istat)

  allocate(xLeft(sizeX),stat=istat)
  allocate(yLeft(sizeY),stat=istat)
  allocate(zLeft(sizeZ),stat=istat)

  allocate(xRight(sizeX),stat=istat)
  allocate(yRight(sizeY),stat=istat)
  allocate(zRight(sizeZ),stat=istat)

  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCenter,sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCenter,sizeZ)

  call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE,gcell,xLeft,sizeX)
  call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE,gcell,yLeft,sizeY)
  call Grid_getCellCoords(KAXIS,blockID,LEFT_EDGE,gcell,zLeft,sizeZ)

  call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,gcell,xRight,sizeX)
  call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,gcell,yRight,sizeY)
  call Grid_getCellCoords(KAXIS,blockID,RIGHT_EDGE,gcell,zRight,sizeZ)

  ! Assert consistency between simulation geometry and grid geometry
  ! Certain combinations are not implemented
  if (NDIM==1) then
     if (sim_geo==PLANAR_SEDOV) then
        if (.not. (geo==CARTESIAN)) call Driver_abortFlash("Mismatch between sim_geo and geoemetry in 1D")
     else if (sim_geo==CYLINDRICAL_SEDOV) then
        if (.not. (geo==CYLINDRICAL)) call Driver_abortFlash("Mismatch between sim_geo and geoemetry in 1D")   
     else if (sim_geo==SPHERICAL_SEDOV) then
        if (.not. (geo==SPHERICAL)) call Driver_abortFlash("Mismatch between sim_geo and geoemetry in 1D")
     endif 
  endif
  if (NDIM==2) then
     if (sim_geo==PLANAR_SEDOV) then
        call Driver_abortFlash("sim_geo=PLANAR_SEDOV not available in 2D")
     else if (sim_geo==CYLINDRICAL_SEDOV) then
        if (.not. ((geo==CARTESIAN) .or. (geo==POLAR))) &
             call Driver_abortFlash("Must use CARTESIAN or POLAR geometry for sim_geo=CYLINDRICAL_SEDOV in 2D")   
     else if (sim_geo==SPHERICAL_SEDOV) then
        if (.not. ((geo==CYLINDRICAL) .or. (geo==SPHERICAL))) &
             call Driver_abortFlash("Must use CYLINDRICAL or SPHERICAL geometry for sim_geo=SPHERICAL_SEDOV in 2D")
     endif 
  endif
  if (NDIM==3) then
     if (.not. (sim_geo==SPHERICAL_SEDOV)) call Driver_abortFlash("Must use sim_geo=SPHERICAL_SEDOV in 3D")
     if (.not. (geo==CARTESIAN)) call Driver_abortFlash("Must use CARTESIAN geometry in 3D") 
  endif

  ! Must center energy source at origin for polar and spherical
  if ((geo==POLAR) .or. (geo==SPHERICAL)) then
     if ((sim_xcIn>1.0e-7) .or. (sim_ycIn>1.0e-7)) call Driver_abortFlash("Energy source must be at origin for this geometry")
  endif

  ! If using total energy to define source then compute source pressure
  if (sim_useE) then      
     if (sim_geo==PLANAR_SEDOV) Vol = 2.0*sim_RIn
     if (sim_geo==CYLINDRICAL_SEDOV) Vol = PI*sim_RIn**2    
     if (sim_geo==SPHERICAL_SEDOV) Vol = (4.0*PI*sim_RIn**3) / 3.0  

     sim_pIn = (sim_gammaIn - 1.0)*(sim_EIn / Vol)
  endif

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     axis(KAXIS)=k
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS)=j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS)=i

           ! Compute distance to origin for center of cell and each cell corner

           r = (xCenter(i) - sim_xcIn)**2
           r_bsw = (xLeft(i) - sim_xcIn)**2
           r_bse = (xRight(i) - sim_xcIn)**2
           if (NDIM>1) then
              r = r  + (yCenter(j) - sim_ycIn)**2

              r_bne = r_bse
              r_bnw = r_bsw

              r_bne = r_bne + (yRight(j) - sim_ycIn)**2          
              r_bnw = r_bnw  + (yRight(j) - sim_ycIn)**2  
              r_bsw = r_bsw + (yLeft(j) - sim_ycIn)**2  
              r_bse = r_bse + (yLeft(j) - sim_ycIn)**2  
           endif 
           if (NDIM>2) then 
              r = r + (zCenter(k) - sim_zcIn)**2     

              r_tne = r_bne         
              r_tnw = r_bnw
              r_tsw = r_bsw
              r_tse = r_bse  

              r_bne = r_bne + (zLeft(k) - sim_zcIn)**2            
              r_bnw = r_bnw + (zLeft(k) - sim_zcIn)**2  
              r_bsw = r_bsw + (zLeft(k) - sim_zcIn)**2  
              r_bse = r_bse + (zLeft(k) - sim_zcIn)**2   

              r_tne = r_tne + (zRight(k) - sim_zcIn)**2           
              r_tnw = r_tnw + (zRight(k) - sim_zcIn)**2  
              r_tsw = r_tsw + (zRight(k) - sim_zcIn)**2  
              r_tse = r_tse + (zRight(k) - sim_zcIn)**2  
           endif

           r = sqrt(r)
           r_bsw = sqrt(r_bsw)
           r_bse = sqrt(r_bse)
           if (NDIM>1) then
              r_bne = sqrt(r_bne)
              r_bnw = sqrt(r_bnw)
           endif
           if (NDIM>2) then
              r_tne = sqrt(r_tne)
              r_tnw = sqrt(r_tnw)
              r_tsw = sqrt(r_tsw)
              r_tse = sqrt(r_tse)
           endif

           if ((geo==POLAR) .or. (geo==SPHERICAL)) then
              r_bne = xRight(i)
              r_bnw = xLeft(i)
              r_bsw = xLeft(i)
              r_bse = xRight(i)
           endif


           ! 1D Calculation

           if (NDIM==1) then
              if (r_bse <= sim_rIn) then
              ! Inside energy source
                 mfrac1 = 1.0
                 mfrac2 = 0.0
                 mfrac3 = 0.0
                 rho = sim_rhoIn
                 p   = sim_pIn
              else if (r_bsw > sim_rIn) then
              ! Outside energy source
                 if (sim_ibound .and. (xCenter(i) > sim_h1)) then
                 ! In 2nd atmosphere
                    mfrac1 = 0.0
                    mfrac2 = 0.0
                    mfrac3 = 1.0
                    call calc_atmos(sim_atmos2,sim_rho2,sim_p2,sim_gamma2,sim_sh2,xCenter(i)-sim_h1,rho,p)
                 else
                 ! In 1st atmosphere
                    mfrac1 = 0.0
                    mfrac2 = 1.0
                    mfrac3 = 0.0
                    call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,xCenter(i),rho,p)
                 endif
              else 
              ! Subsycle computation, boundary of energy source inside the current cell
              ! Density computed using mass fractions, pressure computed using volume averaging
                 dx = (xRight(i) - xLeft(i)) / 10.0
                 pmassIn  = 0.0
                 pmassOut = 0.0
                 p_sum    = 0.0

                 do ii=1,10
                    xx = (ii-0.5)*dx + xLeft(i)
                    rr = sqrt((xx-sim_xcIn)**2) 
                          
                    if (geo==CARTESIAN) dVol = dx

                    dr = dx
                    r2 = ii*dr + xLeft(i)
                    r1 = (ii-1)*dr + xLeft(i)

                    if (geo==CYLINDRICAL) dVol = PI*(r2**2 - r1**2)
                    if (geo==SPHERICAL) dVol = 4.0*PI*(r2**3 - r1**3)/3.0
                         
                    if (rr < sim_rIn) then                            
                       rho_temp = sim_rhoIn
                       p_temp   = sim_pIn
                       p_sum = p_sum + p_temp

                       pmassIn = pmassIn + rho_temp*dVol
                    else 
                       call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,rr,rho_temp,p_temp)
                       p_sum = p_sum + p_temp
                       pmassOut = pmassOut + rho_temp*dVol
                    endif
                    
                 enddo

                 r1 = xLeft(i)
                 r2 = xRight(i)
                 if (geo==CARTESIAN) dVol = r2 - r1
                 if (geo==CYLINDRICAL) dVol = PI*(r2**2 - r1**2)
                 if (geo==SPHERICAL) dVol = 4.0*PI*(r2**3 - r1**3)/3.0    

                 rho = (pmassIn + pmassOut) / dVol
                 p   = p_sum / 10.0

                 mfrac1 = pmassIn / (pmassIn + pmassOut)
                 mfrac2 = pmassOut / (pmassIn + pmassOut)
                 mfrac3 = 0.0
              endif                      
           endif



           if (NDIM==2) then
              if ((r_bne <= sim_rIn) .and. (r_bnw <= sim_rIn) .and. (r_bsw <= sim_rIn) .and. (r_bse <= sim_rIn)) then
              ! Inside energy source
                 mfrac1 = 1.0
                 mfrac2 = 0.0
                 mfrac3 = 0.0
                 rho = sim_rhoIn
                 p   = sim_pIn
              else
                 if ((r_bne > sim_rIn) .and. (r_bnw > sim_rIn) .and. (r_bsw > sim_rIn) .and. (r_bse > sim_rIn)) then
                    if ((geo==POLAR) .or. (geo==SPHERICAL)) then
                       ! Outside energy source
                       if ((sim_ibound) .and. (xCenter(i) > sim_h1)) then
                       ! In 2nd atmosphere
                          mfrac1 = 0.0
                          mfrac2 = 0.0
                          mfrac3 = 1.0    
                          call calc_atmos(sim_atmos2,sim_rho2,sim_p2,sim_gamma2,sim_sh2,xCenter(i)-sim_h1,rho,p)   
                       else 
                       ! In 1st atmosphere   
                          mfrac1 = 0.0
                          mfrac2 = 1.0
                          mfrac3 = 0.0
                          call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,xCenter(i),rho,p)
                       endif
                    else
                       ! Outside energy source
                       if ((sim_ibound) .and. (yCenter(j) > sim_h1)) then
                       ! In 2nd atmosphere
                          mfrac1 = 0.0
                          mfrac2 = 0.0
                          mfrac3 = 1.0    
                          call calc_atmos(sim_atmos2,sim_rho2,sim_p2,sim_gamma2,sim_sh2,yCenter(j)-sim_h1,rho,p)   
                       else 
                       ! In 1st atmosphere   
                          mfrac1 = 0.0
                          mfrac2 = 1.0
                          mfrac3 = 0.0
                          call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,yCenter(j),rho,p)
                       endif
                    endif
                 else 
                 ! Subsycle computation, boundary of energy source inside the current cell
                 ! Density computed using mass fractions, pressure computed using volume averaging

                    dx = (xRight(i) - xLeft(i)) / 10.0
                    dy = (yRight(j) - yLeft(j)) / 10.0

                    pmassIn  = 0.0
                    pmassOut = 0.0
                    p_sum    = 0.0

                    do ii=1,10
                       do jj=1,10
                          xx = (ii-0.5)*dx + xLeft(i)
                          yy = (jj-0.5)*dy + yLeft(j)
                          rr = sqrt((xx-sim_xcIn)**2 + (yy-sim_ycIn)**2)
                          if ((geo==POLAR) .or. (geo==SPHERICAL)) rr = xx
                                 
                          if (geo==CYLINDRICAL) then
                             r1 = (ii-1)*dx + xLeft(i)
                             r2 = ii*dx + xLeft(i)
                             dVol = PI*(r2**2 - r1**2)*dy
                          else if (geo==POLAR) then
                             r1 = (ii-1)*dx + xLeft(i)
                             r2 = ii*dx + xLeft(i)
                             dVol = 0.5*(r2**2 - r1**2)*dy       
                          else if (geo==SPHERICAL) then 
                             r1 = (ii-1)*dx + xLeft(i)
                             r2 = ii*dx + xLeft(i)
                             theta1 = (jj-1)*dy + yLeft(j)
                             theta2 = jj*dy + yLeft(j)   
                             dVol = 2.0*PI*(r2**3 - r1**3)*(cos(theta1) - cos(theta2))/3.0                                             
                          else
                             dVol = dx*dy
                          endif          
                         
                          if (rr < sim_rIn) then                            
                             rho_temp = sim_rhoIn
                             p_temp   = sim_pIn
                             p_sum = p_sum + p_temp

                             pmassIn = pmassIn + rho_temp*dVol
                          else 
                             if ((geo==POLAR) .or. (geo==SPHERICAL)) then
                                call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,rr,rho_temp,p_temp)
                             else 
                                call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,yy,rho_temp,p_temp)
                             endif 

                             p_sum = p_sum + p_temp
                             pmassOut = pmassOut + rho_temp*dVol
                          endif
                       enddo
                    enddo

                    dx = xRight(i) - xLeft(i)
                    dy = yRight(j) - yLeft(j)

                    if (geo==CYLINDRICAL) then
                       r1 = xLeft(i)
                       r2 = xRight(i)
                       dVol = PI*(r2**2 - r1**2)*dy
                    else if (geo==POLAR) then
                       r1 = xLeft(i)
                       r2 = yRight(i)
                       dVol = 0.5*(r2**2 - r1**2)*dy      
                    else if (geo==SPHERICAL) then
                       r1 = xLeft(i)
                       r2 = xRight(i)
                       theta1 = yLeft(j)
                       theta2 = yRight(j)
                       dVol = 2.0*PI*(r2**3 - r1**3)*(cos(theta1) - cos(theta2))/3.0  
                    else
                       dVol = dx*dy
                    endif          

                    rho = (pmassIn + pmassOut) / dVol
                    p   = p_sum / 100.0

                    mfrac1 = pmassIn / (pmassIn + pmassOut)
                    mfrac2 = pmassOut / (pmassIn + pmassOut)
                    mfrac3 = 0.0
                 endif
              endif
           endif



           if (NDIM==3) then
              if ((r_bne <= sim_rIn) .and. (r_bnw <= sim_rIn) .and. (r_bsw <= sim_rIn) .and. (r_bse <= sim_rIn) .and. &
                   (r_tne <= sim_rIn) .and. (r_tnw <= sim_rIn) .and. (r_tsw <= sim_rIn) .and. (r_tse <= sim_rIn)) then
              ! Inside energy source
                 mfrac1 = 1.0
                 mfrac2 = 0.0
                 mfrac3 = 0.0
                 rho = sim_rhoIn
                 p   = sim_pIn
              else
                 if ((r_bne > sim_rIn) .and. (r_bnw > sim_rIn) .and. (r_bsw > sim_rIn) .and. (r_bse > sim_rIn) .and. &
                      (r_tne > sim_rIn) .and. (r_tnw > sim_rIn) .and. (r_tsw > sim_rIn) .and. (r_tse > sim_rIn)) then
                    ! Outside energy source
                    if ((sim_ibound) .and. (zCenter(j) > sim_h1)) then
                    ! In 2nd atmosphere
                       mfrac1 = 0.0
                       mfrac2 = 0.0
                       mfrac3 = 1.0    
                       call calc_atmos(sim_atmos2,sim_rho2,sim_p2,sim_gamma2,sim_sh2,zCenter(j)-sim_h1,rho,p)   
                    else 
                    ! In 1st atmosphere   
                       mfrac1 = 0.0
                       mfrac2 = 1.0
                       mfrac3 = 0.0
                       call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,zCenter(j),rho,p)
                    endif
                 else
                    
! write (*,*) 'Subway'
! write (*,*) 'xyz=', xCenter(i),yCenter(j),zCenter(k)

                 ! Subsycle computation, boundary of energy source inside the current cell
                 ! Density computed using mass fractions, pressure computed using volume averaging

                    dx = (xRight(i) - xLeft(i)) / 10.0
                    dy = (yRight(j) - yLeft(j)) / 10.0
                    dz = (zRight(k) - zLeft(k)) / 10.0                  

                    pmassIn  = 0.0
                    pmassOut = 0.0
                    p_sum    = 0.0

                    do ii=1,2
                       do jj=1,2
                          do kk=1,2
                             xx = (ii-0.5)*dx + xLeft(i)
                             yy = (jj-0.5)*dy + yLeft(j)
                             zz = (kk-0.5)*dz + ZLeft(k)

                             rr = sqrt((xx-sim_xcIn)**2 + (yy-sim_ycIn)**2 + (zz-sim_zcIn)**2)

                             dVol = dx*dy*dz
                     
                             if (rr < sim_rIn) then                            
                                rho_temp = sim_rhoIn
                                p_temp   = sim_pIn
                                p_sum = p_sum + p_temp

                                pmassIn = pmassIn + rho_temp*dVol
                             else 
                                call calc_atmos(sim_atmos1,sim_rho1,sim_p1,sim_gamma1,sim_sh1,zz,rho_temp,p_temp)
                                p_sum = p_sum + p_temp
                                pmassOut = pmassOut + rho_temp*dVol
                             endif
                          enddo
                       enddo
                    enddo

                    dx = xRight(i) - xLeft(i)
                    dy = yRight(j) - yLeft(j)
                    dz = zRight(k) - zLeft(k)

                    dVol = dx*dy*dz

                    rho = (pmassIn + pmassOut) / dVol
                    p   = p_sum / 1000.0

                    mfrac1 = pmassIn / (pmassIn + pmassOut)
                    mfrac2 = pmassOut / (pmassIn + pmassOut)
                    mfrac3 = 0.0
                 endif
              endif
           endif

           call Grid_putPointData(blockID,CENTER,VELX_VAR,EXTERIOR,&
                                  axis,0.0)
           call Grid_putPointData(blockID,CENTER,VELY_VAR,EXTERIOR,&
                                  axis,0.0)
           call Grid_putPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,&
                                  axis,0.0)
           call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,&
                                  axis,rho)
           call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,&
                                  axis,p)
!           call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,&
!                                  axis,eint)
           call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN,EXTERIOR,&
                                  axis,mfrac1)
           call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+1,EXTERIOR,&
                                  axis,mfrac2)
           call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+2,EXTERIOR,&
                                  axis,mfrac3)
        enddo
     enddo
  enddo





! Thermodynamic calculations
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     axis(KAXIS)=k
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS)=j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS)=i

           call Grid_getPointData(blockID,CENTER,DENS_VAR,EXTERIOR,&
                                  axis,rho)
           call Grid_getPointData(blockID,CENTER,PRES_VAR,EXTERIOR,&
                                  axis,p)
 !          call Grid_getPointData(blockID,CENTER,EINT_VAR,EXTERIOR,&
 !                                 axis,eint)
           call Grid_getPointData(blockID,CENTER,SPECIES_BEGIN,EXTERIOR,&
                                  axis,mfrac1)
           call Grid_getPointData(blockID,CENTER,SPECIES_BEGIN+1,EXTERIOR,&
                                  axis,mfrac2)
           call Grid_getPointData(blockID,CENTER,SPECIES_BEGIN+2,EXTERIOR,&
                                  axis,mfrac3)

           eosData(EOS_DENS) = rho
           eosData(EOS_PRES) = p
!          eosData(EOS_EINT) = eint

           massFrac(1) = mfrac1
           massFrac(2) = mfrac2
           massFrac(3) = mfrac3 

!           p = mfrac1*sim_pIn + mfrac2*sim_p1 + mfrac3*sim_p2
           eosData(EOS_PRES) = p

           call Eos(MODE_DENS_PRES,1,eosData,massFrac,mask)

           eint = eosData(EOS_EINT)
           temp = eosData(EOS_TEMP)

           gamc = mfrac1*sim_gammaIn + mfrac2*sim_gamma1 + mfrac3*sim_gamma2
           game = gamc

           call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,&
                                  axis,p)
           call Grid_putPointData(blockID,CENTER,GAMC_VAR,EXTERIOR,&
                                  axis,gamc)
           call Grid_putPointData(blockID,CENTER,GAME_VAR,EXTERIOR,&
                                  axis,game)
           call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,&
                                  axis,eint)
           call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,&
                                  axis,eint)
           call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,&
                                  axis,temp)
        enddo
     enddo
  enddo
  !     
  
  return
end subroutine Simulation_initBlock
