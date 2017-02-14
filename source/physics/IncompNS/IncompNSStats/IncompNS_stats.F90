!!****if* source/physics/IncompNS/IncompNSStats/IncompNS_stats
!!
!! NAME
!!
!!  IncompNS_stats
!!
!! SYNOPSIS
!!
!!  call IncompNS_stats()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***



subroutine IncompNS_stats()

  use Driver_interface, only : Driver_getSimTime
  use IncompNS_data, only : ins_meshMe,ins_nstep, ins_useIncompNS
  use ins_statsData, only : ins_statsN,ins_statsIntervalStep,ins_statsStartTime
  use ins_statsInterface, only : ins_statsVelpTimeAvg, ins_statsRestressesTimeavg 

  use Grid_interface, ONLY : Grid_computeVarMean

  implicit none
#include "constants.h"
#include "Flash.h"

  logical :: stats_flg
  real    :: simTime

  real :: w,w2

  if (.NOT. ins_useIncompNS) RETURN

  stats_flg = (MOD(ins_nstep,ins_statsIntervalStep) .eq. 0)

  if (stats_flg) then

    call Driver_getSimTime(simTime)

    if (simTime .ge. ins_statsStartTime) then

    ! Velocities Time average
    call ins_statsVelpTimeAvg(ins_statsN)

    ! Reynolds Stresses Time Average
    call ins_statsRestressesTimeavg(ins_statsN)

    ! etc.

    ! Update ensemble number:
    ins_statsN = ins_statsN + 1

    ! Monitor w vel:
    call Grid_computeVarMean(DUST_VAR,w)  ! <w> vel that has been interpolated to cell center in ins_statsVelpTimeAvg
    call Grid_computeVarMean(WWAV_VAR,w2) ! <w*w> + <w>*<w>   

    if (ins_meshMe .eq. MASTER_PE) then
       write(*,*) "Time averaged statistics Computed, ensemble num=",ins_statsN
       write(*,*) "Time+domain averaged <w^2>=",w2-w*w,", where mean velocity <w>=",w
    endif

    else

      if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) ' '
        write(*,*) &
        'Statistics delayed, simTime=',simTime,' less than Stats start time=',ins_statsStartTime
      endif

    endif

  endif

  return

end subroutine

