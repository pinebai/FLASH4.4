!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( integer, intent(IN)    :: blockCount, 
!!               integer(:), intent(IN) :: blockList, 
!!               real, intent(IN)       ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
!!   dt  --       passed to the internal bn_burner module  
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useBurnTable -- Boolean, False.  Controls the generation of reaction rates.
!!                TRUE interpolates from a stored table; FALSE generates them
!!                analytically.
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  algebra -- Integer, 1, [1,2].  Controls choice of linear algebra package used
!!                for matrix solution.  1=Ma28 sparse package, 2=Gift hardwired package.
!!  odeStepper -- Integer, 1, [1,2].  Controls time integration routines.
!!                1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!  nuclearTempMin/Max -- Real, 1.1E+8/1.0E+12.  Minimum and maximum temperature
!!                ranges where burning can occur
!!  nuclearDensMin/Max -- Real, 1.0E-10/1.0E+14.  Minimum and maximum density range
!!                where burning can occur.
!!  nuclearNI56Max -- Real, 1.0.  Maximum mass fraction of nickel where burning
!!                can occur.
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.See Burn_computeDt for details.              
!!
!! NOTES
!!
!!  The burning unit adds a new mesh variable ENUC_VAR which is the nuclear energy 
!!             generation rate
!!
!!***

!!REORDER(4): solnData

subroutine Burn (  blockCount, blockList, dt  )    

  use Burn_data, ONLY:  bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_smallx, bn_useBurn, bn_meshMe
  use bn_interface, ONLY :  bn_mapNetworkToSpecies, bn_burner   


  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Hydro_interface, ONLY:  Hydro_detectShock

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  !args
  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  ! locals
  integer                 :: i, j, k, n, specieMap
  integer                 :: blockID, thisBlock

  real, dimension(NSPECIES) :: xIn, xOut
  real                    :: sdot
  real                    :: tmp, rho, ei, ek
  logical                 :: burnedZone

  integer, parameter      :: NNVARS = 11
  integer, parameter      :: nin  = NSPECIES+NNVARS 
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  logical okBurnTemp, okBurnDens, okBurnShock, okBurnNickel  ! determines when Burning actually occurs.  
  logical :: getGuardCells = .true.
  !!     Sorry, can no longer do this trick with interfaces defined.  
  !!       Need to allocated them the size of xSizeCoord etc.
  !!       real, dimension(MAXCELLS) :: xCoord, yCoord, zCoord
  real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord
  integer                                 :: xSizeCoord, ySizeCoord, zSizeCoord

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: shock
#endif

  real, pointer, dimension(:,:,:,:)            :: solnData

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

  !---------------------------------------------------------------------------------

  ! start the timer ticking
  call Timers_start("burn")

  ! make sure that guardcells are up to date
  if (.NOT. bn_useShockBurn) then


     call Grid_fillGuardCells(CENTER, ALLDIR)

  endif

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount

     blockID = blockList(thisBlock)
     burnedZone = .FALSE.

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     xSizeCoord = blkLimitsGC(HIGH,IAXIS)
     ySizeCoord = blkLimitsGC(HIGH,JAXIS)
     zSizeCoord = blkLimitsGC(HIGH,KAXIS)
     !! allocate space for dimensions
     allocate(xCoord(xSizeCoord))
     allocate(yCoord(ySizeCoord))
     allocate(zCoord(zSizeCoord))

     call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockID,solnData)

     if (.NOT. bn_useShockBurn) then
        call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
             xCoord,yCoord,zCoord)
     else
        shock(:,:,:) = 0
     endif


     ! now guaranteed that tmp, rho, etc. exist
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              okBurnTemp = .FALSE.
              okBurnDens = .FALSE.
              okBurnShock = .FALSE.
              okBurnNickel = .FALSE.

              tmp  = solnData(TEMP_VAR,i,j,k)
              rho  = solnData(DENS_VAR,i,j,k)

              sdot = 0.0e0

  ! Split up this monster IF statement so debuggers can tell when it is passing/failing
  !            if ( (tmp >= bn_nuclearTempMin .AND. tmp <= bn_nuclearTempMax) .AND. & 
  !                 (rho >= bn_nuclearDensMin .AND. rho <= bn_nuclearDensMax) .AND. & 
  !                 (shock(i,j,k) == 0.0 .OR.  & 
  !                 (shock(i,j,k) == 1.0 .AND. bn_useShockBurn))) then
              okBurnTemp = (tmp >= bn_nuclearTempMin .AND. tmp <= bn_nuclearTempMax)
              okBurnDens = (rho >= bn_nuclearDensMin .AND. rho <= bn_nuclearDensMax) 
              okBurnShock = (shock(i,j,k) == 0.0 .OR. (shock(i,j,k) == 1.0 .AND. bn_useShockBurn))
              if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

  !  Again, split up this if statement into something LBR can understand and test
  !                 if ( (NI56_SPEC <= 0) .OR.  & 
  !                      (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max) ) then
                 if (NI56_SPEC /= NONEXISTENT) then            
                    okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                 else    ! nickel is not even a species in this simulation, so we'll always burn
                    okBurnNickel = .TRUE.
                 endif

                 if (okBurnNickel) then

                    burnedZone = .TRUE.

                    ! Map the solution data into the order required by bn_burner
                    do n = 1, NSPECIES
                       call bn_mapNetworkToSpecies(n,specieMap)
                       xin(n) = solnData(specieMap,i,j,k)
                    end do

                     

                    ! Do the actual burn
                    call bn_burner(dt, tmp, rho, xIn, xOut, sdot)

                    !  Map it back NOTE someday make a nicer interface....
                    do n = 1, NSPECIES
                       call bn_mapNetworkToSpecies(n,specieMap)
                       solnData(specieMap,i,j,k) = xOut(n)
                    end do

   !  NOTE should probably do something here with eintSwitch for consistency
   !  LBR will settle for simply using internal energy!
                    ! kinetic energy
                    ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  & 
                         solnData(VELY_VAR,i,j,k)**2 +  & 
                         solnData(VELZ_VAR,i,j,k)**2)

                    ! internal energy, add on nuclear rate*timestep
                    ei = solnData(ENER_VAR,i,j,k) - ek
                    ei = ei + dt*sdot
                     
 
#ifdef EINT_VAR
                    solnData(EINT_VAR,i,j,k) = ei
#endif
                    solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                    solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + dt*sdot
#endif
                    solnData(ENUC_VAR,i,j,k) = sdot

                 endif
              endif

           enddo
        enddo
     enddo

     ! we've altered the EI, let's equilabrate
     if (burnedZone) then

#ifdef FLASH_UHD_3T
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID) ! modified for 3T
#else 
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
#endif

        !Notes from flash2 -- new eos below seems to work for almost every setup
        !except hse_isothermal_fuel+ash. no idea why. In the
        !meantime we'll use old verstion
        !          call eos3d_new(solnData,0,eos_order,               &
        !                      nxb,    nyb,    nzb,               &
        !                      k2d,    k3d,    nguard,            &
        !                      maxcells, NSPECIES, nin,             &
        !                      eos_params)
     end if


     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)

  end do

  call Timers_stop("burn")

  return
end subroutine Burn
