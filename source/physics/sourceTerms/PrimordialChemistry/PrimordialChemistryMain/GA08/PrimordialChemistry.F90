!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/PrimordialChemistry
!!
!! NAME
!!
!! PrimordialChemistry
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry( integer         :: blockCount,
!!                            integer(:)      :: blockList,
!!                            real            :: dt)
!!
!! DESCRIPTION
!!
!! Apply chemistry to all blocks in specified list
!!
!! ARGUMENTS
!!
!! blockCount -- dimension of blockList
!! blockList  -- array of blocks which should receive chemistry
!! dt   --       passed to the internal pchem_burner module
!!
!! PARAMETERS
!!
!! usePrimordialChemistry  -- Boolean, True. Turns on chemistry module
!! useShockChem  -- ?? Nothing for now, may need later
!! pchem_algebra       -- Integer, 2, should control choice of linear pchem_algebra 
!!                  packages, 1-MA28 for sparse matricies, 2-GIFT for
!!                  non sparse. HAVE TO USE 2 for CHEMISTRY, 1 is not supported
!! pchem_odeStepper    -- Integer, 1, [1,2]. Controls time integration routines
!!                  1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!
!! NOTES
!!
!!
!!***

subroutine PrimordialChemistry(blockCount, blockList, dt)
   
   use PrimordialChemistry_data
   use pchem_interface
   use Simulation_data

   use Timers_interface, ONLY : Timers_start, Timers_stop
   use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getBlkIndexLimits, &
                              Grid_getCellCoords, Grid_getBlkPtr, &
                              Grid_releaseBlkPtr
   use Eos_interface, ONLY : Eos_wrapped, Eos
   use Hydro_interface, ONLY : Hydro_detectShock

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

   !args
   integer, INTENT(in)                          :: blockCount
   integer, INTENT(in), DIMENSION(blockCount)   :: blockList
   real, intent(IN)                             :: dt

   !locals
   integer                                      :: i, j, k, n, specieMap
   integer                                      :: blockID, thisBlock

   real, dimension(NSPECIES)    :: xIn, xOut
   real                         :: sdot
   real                         :: tmp, rho, ei, ek, temp_back
   logical                      :: chemZone
   integer                      :: sizex, mcheck

   integer, parameter           :: NNVARS = 11 !what is this?
   integer, parameter           :: nin = NSPECIES+NNVARS
   integer, dimension(2,MDIM)   :: blkLimits, blkLimitsGC

   logical okChemTemp, okChemDens, okChemShock, okChemNickel !! Maynot need these
   logical :: getGuardCells = .true.

   real, allocatable, dimension(:)      :: xCoord, yCoord, zCoord
   integer                              :: xSizeCoord, ySizeCoord, zSizeCoord
   
   real                                 :: counts, jcounts
   real                                 :: metal_frac
   real                                 :: eosData(EOS_NUM),xeos(SPECIES_BEGIN:SPECIES_END)


#ifdef FIXEDBLOCKSIZE
      real, dimension(GRID_IHI_GC, GRID_JHI_GC, GRID_KHI_GC) :: shock
#endif



   real, pointer, dimension(:,:,:,:)    :: solnData

   ! -------------------- Check if chemistry is requested in runtime parameter
   if (.not. pchem_usePrimordialChemistry) return

   !--- Off to the races, 

   call Timers_start("chemistry")
   
   if (.not.  pchem_useShockBurn) then
      call Grid_fillGuardCells(pchem_meshMe,CENTER,ALLDIR)
   endif


!   print *, 'RCCASE: ', pchem_rcCase
!   print *, 'CCCASE: ', pchem_ccCase
!   print *, 'J21: ', pchem_j21

   ! loop over list of blocks passed in
   do thisBlock = 1, blockCount
        
      blockID = blockList(thisBlock)
      !chemZone = .FALSE. !Don't think I need this since we want Chem everwhere
      ! get dimensions/limits and coordinates
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      xSizeCoord = blkLimitsGC(HIGH,IAXIS)
      ysizeCoord = blkLimitsGC(HIGH,JAXIS)
      zsizeCoord = blkLimitsGC(HIGH,KAXIS)
      !! Allocate space for dimensions
      allocate(xCoord(xSizeCoord))
      allocate(yCoord(ySizeCoord))
      allocate(zCoord(zSizeCoord))

      call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
      call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
      call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

      ! Get a pointer to solution data
      call Grid_getBlkPtr(blockID,solnData)

!!$      if (.NOT. pchem_useShockBurn) then
!!$         call Hydro_detectShock(solnData,shock,blkLimits,blkLimitsGC, 1, &
!!$                             xCoord,yCoord,zCoord)
!!$      else
!!$     shock(:,:,:) = 0
!!$      endif

      mcheck = 0     !!Assume no metals
#ifdef METL_MSCALAR
      mcheck = 1  !! if metl_mscalar exists we use it for metals
#endif

      ! now guaranteed that tmp, rho, etc. exist
      do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
         do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                  okChemTemp = .FALSE.
                  okChemDens = .FALSE.
                  okChemShock = .FALSE.
  
                 tmp = solnData(TEMP_VAR,i,j,k)
                 rho = solnData(DENS_VAR,i,j,k)
                 if(pchem_mCool .eq. 1 .and. mcheck .eq. 1) then
                   metal_frac = solnData(METL_MSCALAR,i,j,k)/rho
                 else if(pchem_mCool .eq. 0 .and. mcheck .eq. 1) then
                   metal_frac = 0.0
                 else !!catch all. just put mfrac at zero
                   metal_frac = 0.0
                 endif

!               print *, 'CHEMISTRY:   TEMP: ' , tmp, ' DENS: ', rho
!                 tmp= sim_c_temp
!                 rho= sim_c_den

                  sdot = 0.0e0
                  temp_back = 0.0e0
                  ek = 0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                            solnData(VELY_VAR,i,j,k)**2 + &
                            solnData(VELZ_VAR,i,j,k)**2 )
                  ei = solnData(EINT_VAR,i,j,k)

              

  !! The burn unit has a lot of if statements here. But they are checking a lot of
  !! things to see if they should turn burning on in a certain block. I don't think
  !! I have to worry about this, since we want the chemistry to run everywhere.

                  do n = 1, NSPECIES
                     call pchem_mapNetworkToSpecies(n,specieMap)
                     xIn(n) = solnData(specieMap,i,j,k)
                         !  print *, 'xin(',n,')=',xIn(n), 'n:' , n
                  enddo

                  if(tmp .lt. 50.0) then
                     tmp = 51.0
                     xeos(H_SPEC) = xIn(iH)
                     xeos(HP_SPEC) = xIn(iHP)
                     xeos(HM_SPEC) = xIn(iHM)
                     xeos(D_SPEC) = xIn(iD)
                     xeos(DP_SPEC) = xIn(iDP)
                     xeos(DM_SPEC) = xIn(iDM)
                     xeos(HE_SPEC) = xIn(iHE)
                     xeos(HEP_SPEC) = xIn(iHEP)
                     xeos(HEPP_SPEC) = xIn(iHEPP)
                     xeos(H2_SPEC) = xIn(iH2)
                     xeos(H2P_SPEC) = xIn(iH2P)
                     xeos(HD_SPEC) = xIn(iHD)
                     xeos(HDP_SPEC) = xIn(iHDP)
                     xeos(D2_SPEC) = xIn(iD2)
                     xeos(D2P_SPEC) = xIn(iD2P)
                     xeos(ELEC_SPEC) = xIn(iELEC) 
                     eosData(EOS_TEMP) = tmp
                     eosData(EOS_DENS) = rho
                     eosData(EOS_EINT) = ei
                     call Eos(MODE_DENS_TEMP,1,eosData,xeos)
                     tmp = eosData(EOS_TEMP)
                     rho = eosData(EOS_DENS)
                     ei = eosData(EOS_EINT)
                  endif
                  
                
                  !Do the chemistry

                  call pchem_burner(dt,tmp,rho,xIn,xOut,sdot,ei,counts, jcounts, metal_frac,temp_back)


                  !map the species back
                  do n=1,NSPECIES
                     call pchem_mapNetworkToSpecies(n,specieMap)
                     solnData(specieMap,i,j,k) = xOut(n)
                  enddo
        

                  if(temp_back .lt. 50.0) then
                     temp_back = 51.0
                     xeos(H_SPEC) = xOut(iH)
                     xeos(HP_SPEC) = xOut(iHP)
                     xeos(HM_SPEC) = xOut(iHM)
                     xeos(D_SPEC) = xOut(iD)
                     xeos(DP_SPEC) = xOut(iDP)
                     xeos(DM_SPEC) = xOut(iDM)
                     xeos(HE_SPEC) = xOut(iHE)
                     xeos(HEP_SPEC) = xOut(iHEP)
                     xeos(HEPP_SPEC) = xOut(iHEPP)
                     xeos(H2_SPEC) = xOut(iH2)
                     xeos(H2P_SPEC) = xOut(iH2P)
                     xeos(HD_SPEC) = xOut(iHD)
                     xeos(HDP_SPEC) = xOut(iHDP)
                     xeos(D2_SPEC) = xOut(iD2)
                     xeos(D2P_SPEC) = xOut(iD2P)
                     xeos(ELEC_SPEC) = xOut(iELEC) 
                     eosData(EOS_TEMP) = temp_back
                     eosData(EOS_DENS) = rho
                     eosData(EOS_EINT) = ei
                     call Eos(MODE_DENS_TEMP,1,eosData,xeos)
                     temp_back = eosData(EOS_TEMP)
                     ei = eosData(EOS_EINT)
                     rho = eosData(EOS_DENS)
                  endif

                     solnData(TEMP_VAR,i,j,k) = temp_back
                     solnData(DENS_VAR,i,j,k) = rho
                     solnData(EINT_VAR,i,j,k) = ei
                     solnData(CHDT_VAR,i,j,k) = counts
                     solnData(JTDT_VAR,i,j,k) = jcounts
                     solnData(ENER_VAR,i,j,k) = ei+ek

           enddo
        enddo
     enddo


  !we've altered the EI, let's equilabrate
!   call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
   call Eos_wrapped(MODE_DENS_EI,blkLimitsGC,blockID)
!    call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
!   call Eos_wrapped(MODE_DENS_PRES,blkLimits,blockID) 
 
    call Grid_releaseBlkPtr(blockID,solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)
  enddo

  call Timers_stop("chemistry")
  
  return

end subroutine PrimordialChemistry
