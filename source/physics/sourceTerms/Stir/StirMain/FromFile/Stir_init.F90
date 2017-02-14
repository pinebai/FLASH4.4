!!****if* source/physics/sourceTerms/Stir/StirMain/FromFile/Stir_init
!!
!! NAME
!!  Stir_init
!!
!! SYNOPSIS
!!  call Stir_init(logical(in) :: restart)
!!
!! DESCRIPTION
!!  Read forcing pattern from file and apply the stirring operator
!!  on the list of blocks provided as input
!!
!! ARGUMENTS
!!   restart - restarting from checkpoint?
!!
!! PARAMETERS
!!   These are the runtime parameters used in the Stir unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.
!!
!!    useStir        [BOOLEAN]
!!        Switch to turn stirring on or off at runtime.
!!    st_infilename  [CHARACTER]
!!        file containing the stirring time sequence
!!    st_computeDt   [BOOLEAN]
!!        whether to restrict timestep based on stirring
!!
!! AUTHOR
!!  Christoph Federrath, 2008
!!
!!***

subroutine Stir_init(restart)

  use Stir_data
  use Driver_interface,            ONLY : Driver_getSimTime, Driver_abortFlash, &
                                          Driver_getComm, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface,              ONLY : Grid_getGlobalIndexLimits
  use IO_interface,                ONLY : IO_getScalar, IO_setScalar
  implicit none

#include "constants.h"
#include "Flash.h"

  logical, intent(in) :: restart

  real                :: time, timeinfile
  logical, parameter  :: Debug = .false.

  call RuntimeParameters_get('useStir', st_useStir)
  call RuntimeParameters_get('st_infilename', st_infilename)
  call RuntimeParameters_get('st_computeDt', st_computeDt)

  if (.not.st_useStir) write(*,'(A)') 'WARNING:  You have included the StirFromFile unit but useStir=.false.'

  call Driver_getMype(GLOBAL_COMM, st_globalMe)
  call Driver_getMype(MESH_COMM,st_meshMe)
  call Driver_getComm(MESH_COMM,st_meshComm)

  ! this call sets st_dtDpdateAccel and reads general information from file (see below)
  call Driver_getSimTime(time)
  call st_readStirringDataFromFile(st_infilename, time, timeinfile)

  if (st_globalMe.eq.MASTER_PE) then
     write (*,'(A,I4,A)') 'read ',st_nmodes,' modes for stirring from file: ', trim(st_infilename)
     if (st_spectform == 0) write (*,'(A,I2,A)') ' spectral form        = ', st_spectform, ' (Band)'
     if (st_spectform == 1) write (*,'(A,I2,A)') ' spectral form        = ', st_spectform, ' (Paraboloid)'
     write (*,'(A,ES10.3)') ' solenoidal weight    = ', st_solweight
     write (*,'(A,ES10.3)') ' st_solweightnorm     = ', st_solweightnorm
     write (*,'(A,ES10.3)') ' stirring energy      = ', st_energy
     write (*,'(A,ES10.3)') ' autocorrelation time = ', st_decay
     write (*,'(A,ES10.3)') ' minimum wavenumber   = ', st_stirmin
     write (*,'(A,ES10.3)') ' maximum wavenumber   = ', st_stirmax
  endif

  return

end subroutine Stir_init
