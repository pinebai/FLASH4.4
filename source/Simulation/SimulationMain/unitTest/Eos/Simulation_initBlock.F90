!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (in) :: blockId, 
!!                   
!!
!!
!!
!! DESCRIPTION
!!   initializes a 3-d domain by setting with a density 
!!   gradient in x, temperature gradient in y, and composition gradient
!!   in z. An Eos call is made with mode = MODE_DENSITY_TEMP, and the
!!   resultant pressure and energy are saved. 
!!  
!!
!! 
!! ARGUMENTS
!!
!!  blockId -        the ID of block to put initial condtions into
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_zmin,sim_zmax,sim_smallx,&
                              sim_densMin, sim_tempMin, sim_xnMin, sim_presMin,&
                              sim_densMax, sim_tempMax, sim_xnMax, sim_presMax, &
                              sim_initialMass
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_putPointData
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: blockId
  
  
  !
  
  integer :: i, j, k, n, even, put

  integer ::  npts_x, npts_y, npts_z
  real :: xpos, ypos, zpos
  real :: dx, dy, dz

  real ldensMax, ldensMin, dlogrho
  real ltempMax, ltempMin, dlogT
  real lpresMax, lpresMin, dlogP
  real lxnMax, lxnMin, dlogxn
  real :: dens,temp,pres
  real,dimension(NSPECIES) :: xn

  
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  logical,parameter :: gcell = .true.
  integer :: sizeX, sizeY, sizeZ
  integer,dimension(MDIM) :: pos,globalInd
  

  call Grid_getGlobalIndexLimits(globalInd)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,sizeX)
  if(NDIM > 1) then
     call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCenter,sizeY)
  else
     yCenter=0.0
  end if

  if(NDIM > 2) then
     call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCenter,sizeZ)
  else
     zCenter=0.0
  end if

  npts_x = globalInd(IAXIS)
  npts_y = globalInd(JAXIS)
  npts_z = globalInd(KAXIS)
  
! compute the log interval for dens, temp, and xn
  dlogrho = (log10(sim_densMax) - log10(sim_densMin)) / float(npts_x)
  dlogT   = (log10(sim_tempMax) - log10(sim_tempMin)) / float(npts_y)
  dlogP   = (log10(sim_presMax) - log10(sim_presMin)) / float(npts_y)

  if(NSPECIES > 1) then
     dlogxn  = (log10(sim_xnMax)   - log10(sim_xnMin)  ) / float(npts_z)
  end if
! compute the distance between coordinatess -- ideally xmin=0, xmax=1
  dx = (sim_xmax - sim_xmin) / float(npts_x)
  dy = (sim_ymax - sim_ymin) / float(npts_y)
  dz = (sim_zmax - sim_zmin) / float(npts_z)

! since we will be making the density, etc. equally spaced in log, compute
! the log of the extrema
  ldensMax = log10(sim_densMax)
  ldensMin = log10(sim_densMin)
  
  ltempMax = log10(sim_tempMax)
  ltempMin = log10(sim_tempMin)
  
  lpresMax = log10(sim_presMax)
  lpresMin = log10(sim_presMin)
  
  lxnMax   = log10(sim_xnMax)
  lxnMin   = log10(sim_xnMin)

! now fill the master arrays
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

     zpos = (zCenter(k) - sim_zmin) / dz
     pos(KAXIS)=k
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ypos = (yCenter(j) - sim_ymin) / dy
        pos(JAXIS)=j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           xpos = (xCenter(i) - sim_xmin) / dx
           pos(IAXIS)=i
! compute the density, temperature, and composition
           dens = 10.e0**(ldensMin + xpos*dlogrho)
           temp = 10.e0**(ltempMin + ypos*dlogT)
           pres = 10.e0**(lpresMin + ypos*dlogP)


           call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,pos,dens)
           call Grid_putPointData(blockID,CENTER,CTMP_VAR,EXTERIOR,pos,temp)
           ! Need this for Grid_init call on MOD_DENS_TEMP, long before Eos_unitTest is called
           call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,pos,temp)  
           call Grid_putPointData(blockID,CENTER,CPRS_VAR,EXTERIOR,pos,pres)

#if NSPECIES > 0
! for the composition, divide the material according to the input parameter
!  sim_initialMass = -1 to have a gradient divided between species 1 and NSPECIES
!  sim_initialMass = 0 to divide evenly
!  sim_initialMass = i to put all mass in the ith species

           if(NSPECIES==1) then
              xn = 1.0
           else
              select case (sim_initialMass)
                 case (0)              !! divide evenly
                    do even=1, NSPECIES
                       xn(even) = 1.0 / NSPECIES
                     enddo
                  case (-1)            !! anshu's original
                     xn(:) = sim_smallx
                     xn(1) =  10.e0**(lxnMin + zpos*dlogxn)
                     xn(NSPECIES) = 1.e0 - xn(1)
                 case default          !!  put in ith variable
                     xn(:) = sim_smallx
                     xn(sim_initialMass) = 1.0
              end select

              
           end if
           !! Now put the data in UNK
           do put=1,NSPECIES
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+put-1,&
                   EXTERIOR,pos,xn(put))
           enddo
#endif

        enddo
     enddo
  enddo

  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)

  return
end subroutine Simulation_initBlock




