!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/Stir_init
!!
!! NAME
!!
!!  Stir_init
!!
!! SYNOPSIS
!!
!!  call Stir_init(
!!            logical(in) :: restart)
!!
!! DESCRIPTION
!!  Apply the isothermal cooling and stirring opperator 
!!  on the list of blocks provided as input
!!
!! ARGUMENTS
!!   
!!  
!!   restart -indicates if run is starting from scratch or restarting from 
!!            checkpoint
!!
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in the Stir unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!    st_decay [REAL]
!!        correlation time for driving
!!    st_energy [REAL]
!!        energy input/mode
!!    st_freq [INTEGER]
!!        frequency of stirring
!!    st_seed [INTEGER]
!!        random number generator seed
!!    st_stirmax [REAL]
!!        maximum stirring *wavenumber*
!!    st_stirmin [REAL]
!!        minimum stirring *wavenumber*
!!    st_computeDt {BOOLEAN]
!!        whether to restrict timestep based on stirring
!!    useStir [BOOLEAN]
!!        Switch to turn stirring on or off at runtime.
!!***

subroutine Stir_init(restart)

  use Stir_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits,Grid_getBlkIndexLimits
  use IO_interface, ONLY : IO_getScalar, IO_setScalar
  use ut_randomInterface, ONLY : ut_randomSeed
  implicit none

#include "constants.h"
#include "Flash.h"

   
  logical, intent(in) :: restart


  integer :: globalIndexLimits(MDIM)
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: ikxmin, ikxmax, ikymin, ikymax, ikzmin, ikzmax
  integer :: ikx, iky, ikz
  character(len=MAX_STRING_LENGTH) :: eosModeString
  real    :: kx, ky, kz, k, Lx, Ly, Lz, twopi
  real, save :: imin, imax, jmin, jmax, kmin, kmax
  
  call RuntimeParameters_get('st_decay',st_decay)


  call Driver_getMype(MESH_COMM,st_meshMe)
  call Driver_getComm(MESH_COMM,st_meshComm)
  call RuntimeParameters_get( 'xmin',imin)
  call RuntimeParameters_get( 'xmax',imax)
  call RuntimeParameters_get( 'ymin',jmin)
  call RuntimeParameters_get( 'ymax',jmax)
  call RuntimeParameters_get( 'zmin',kmin)
  call RuntimeParameters_get( 'zmax',kmax)
  call RuntimeParameters_get( 'eosMode',eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, st_eosMode)
  
  call RuntimeParameters_get ('st_freq', st_freq)
  call RuntimeParameters_get( 'st_decay',st_decay)
  call RuntimeParameters_get( 'st_energy',st_energy)
  call RuntimeParameters_get( 'st_stirmin',st_stirmin)
  call RuntimeParameters_get( 'st_stirmax',st_stirmax)
  call RuntimeParameters_get( 'st_seed',st_seed)
  call RuntimeParameters_get( 'st_reproducible',st_reproducible)
  call RuntimeParameters_get( 'st_saveReproducible',st_saveReproducible)
  
  call RuntimeParameters_get('useStir', st_useStir)
  call RuntimeParameters_get('st_computeDt', st_computeDt)


  if (.not. st_useStir) then
     write(6,*)'WARNING:  You have included the Stir unit but have set '
     write(6,*)'   the runtime parameter useStir to FALSE'
     write(6,*)'   No stirring will occur but Stir_init will continue.'
  end if

  !initialize some variables, allocate randseed 
  
  st_OUvar = sqrt(st_energy/st_decay)
 

  !if we are starting from scratch
  if(.not. restart) then

     if(st_reproducible.or.st_saveReproducible)&
          open(unit=st_randomSaveUnit,file="saved_random_numbers")
     call ut_randomSeed (ut_size = st_seedLen) 
     if (.not. allocated (st_randseed) ) then
        allocate (st_randseed (st_seedLen) )
     endif
     
     
     if (st_meshMe .eq. MASTER_PE) print *, 'seed length = ', st_seedLen



     !  everyone scan the available modes, and decide which are within the
     !  range for stirring
     !
     ikxmin = 0
     ikymin = 0
     ikzmin = 0
     

#if 0  
     call Grid_getGlobalIndexLimits(globalIndexLimits)
     
     ikxmax = globalIndexLimits(IAXIS)
     ikymax = globalIndexLimits(JAXIS)
     ikzmax = globalIndexLimits(KAXIS)
#endif

  
 
   
     ikxmax = 8
     ikymax = 8
     ikzmax = 8
     
  
     Lx = (imax - imin)
     Ly = (jmax - jmin)
     Lz = (kmax - kmin)
     twopi = 2.*PI
  
  
     st_nmodes = 0
     
     do ikx = ikxmin, ikxmax
        kx = twopi * ikx / Lx
        
        do iky = ikymin, ikymax
           ky = twopi * iky / Ly
           
           do ikz = ikzmin, ikzmax
              kz = twopi * ikz / Lz
              
              k = sqrt( kx*kx+ky*ky+kz*kz )
              
              if ((k.ge.st_stirmin).and.(k.le.st_stirmax))then
                 
                 if ((st_nmodes + 2**(NDIM-1)) .gt. st_maxmodes) then
                    
                    if (st_meshMe .eq. MASTER_PE) print *,&
                         'st_nmodes = ', st_nmodes, ' maxstirmodes = ',st_maxmodes
                    
                    call Driver_abortFlash('Too many stirring modes')
                 endif
                 
                 st_nmodes = st_nmodes + 1
                 
                 st_mode(1,st_nmodes) = kx
                 st_mode(2,st_nmodes) = ky
                 st_mode(3,st_nmodes) = kz
                 
#if NDIM > 1
                 st_nmodes = st_nmodes + 1
                 
                 st_mode(1,st_nmodes) = kx
                 st_mode(2,st_nmodes) =-ky
                 st_mode(3,st_nmodes) = kz
#endif
                 
#if NDIM > 2
                 st_nmodes = st_nmodes + 1
                 
                 st_mode(1,st_nmodes) = kx
                 st_mode(2,st_nmodes) = ky
                 st_mode(3,st_nmodes) =-kz
                 
                 st_nmodes = st_nmodes + 1
                 
                 st_mode(1,st_nmodes) = kx
                 st_mode(2,st_nmodes) =-ky
                 st_mode(3,st_nmodes) =-kz
#endif

                 call IO_setScalar("nmodes", st_nmodes)
                 
              endif
              
           enddo
        enddo
     enddo
     
     if (st_meshMe .eq. MASTER_PE) print *,&
          ' Initializing ',st_nmodes,' modes for stirring.'


     ! Everyone, using the same seed, initialize the OU noises for the
     ! nmodes*6 components of the phases. Store seed in randseed
     ! afterward.
     
     call st_ounoiseinit(st_nmodes*6, st_seed, st_OUvar, st_OUphases)
     call ut_randomSeed (ut_get = st_randseed)
     
  else !we are restarting from a checkpoint

     call IO_getScalar("nmodes", st_nmodes)

     if(st_reproducible.or.st_saveReproducible)&
          open(unit=st_randomSaveUnit,file="saved_random_numbers")


     if (st_meshMe .eq. MASTER_PE) then
       print *, 'seed length = ', st_seedLen          ! if MASTER_PE
       print *, 'Random seed = ', st_randseed         ! if MASTER_PE
       print *, 'nmodes = ', st_nmodes                ! if MASTER_PE
     endif
     call ut_randomSeed (ut_put = st_randseed)
  end if


  ! Then convert those into actual Fourier phases:
  call st_calcPhases()


  return
end subroutine Stir_init
