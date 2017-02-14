!!****if* source/physics/IncompNS/IncompNSStats/ins_statsInit
!!
!! NAME
!!
!!  ins_statsInit
!!
!! SYNOPSIS
!!
!!  call ins_statsInit()
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



subroutine ins_statsInit()

  use Driver_interface, only : Driver_getSimTime
  use IncompNS_data, only : ins_restart,ins_nstep
  use ins_statsData, only : ins_statsIntervalStep,ins_statsRestart,&
                                 ins_statsN,ins_statsStartTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  real :: simTime

  ! How many steps between statistics:
  call RuntimeParameters_get("ins_statsSteps",ins_statsIntervalStep) 

  ! Should we restart the statistics?
  call RuntimeParameters_get("ins_statsRestart",ins_statsRestart)

  ! Compute initial time before start stats:
  call RuntimeParameters_get("ins_statsStartTime",ins_statsStartTime)

  ins_statsN = 0

  if (ins_restart .and. ins_statsRestart) then
     call Driver_getSimTime(simTime)
     if (simTime .ge. ins_statsStartTime) then

    ! Here estimate ensemble number. This overestimates the number if ins_statsStartTime > 0.
        ins_statsN = ins_nstep/ins_statsIntervalStep  

     end if
  endif 

  return

end subroutine

