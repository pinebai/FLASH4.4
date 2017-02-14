!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_setErrorBars
!!
!!  NAME 
!!
!!   sim_setErrorBars
!!
!!  SYNOPSIS
!!
!!   sim_setErrorBars ()
!!
!!  DESCRIPTION
!!
!!   This routine sets the error bars for the simulation. These were determined
!!   from previous runs.
!!
!!***

subroutine sim_setErrorBars ()

  use Simulation_data,             ONLY : sim_focalPercentError1,     &
                                          sim_focalPercentError12,    &
                                          sim_powerPercentError1,     &
                                          sim_powerPercentError12,    &
                                          sim_rayFexitPercentError1,  &
                                          sim_rayFexitPercentError12, &
                                          sim_rayPexitPercentError1,  &
                                          sim_rayPexitPercentError12, &
                                          sim_refinementLevel,        &
                                          sim_symmetryTolerance

  implicit none
!
!
!     ...Set the symmetry tolerance, which checks, if the results are expected
!        due to symmetry.
!
!
  sim_symmetryTolerance = 1.e-8
!
!
!     ...Set the upper bound of the focal point percentage errors. Since these
!        depend on the refinement level, each level has its own upper bound.
!        The ones with the extension '1' apply to those rays having one coordinate
!        equal to the tube's center coordinate and which will feel the deflecting
!        force along only one coordinate. The ones with the extension '12'
!        apply to those rays which will feel the deflecting force along two
!        coordinates.
!
!
  sim_focalPercentError1  (1)  = 0.99
  sim_focalPercentError1  (2)  = 3.28
  sim_focalPercentError1  (3)  = 0.56
  sim_focalPercentError1  (4)  = 0.21
  sim_focalPercentError1  (5)  = 0.01
  sim_focalPercentError1  (6)  = 0.06
  sim_focalPercentError1  (7)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError1  (8)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError1  (9)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError1  (10) = 0.01    ! not actual value, temporary placeholder

  sim_focalPercentError12 (1)  = 7.04
  sim_focalPercentError12 (2)  = 4.56
  sim_focalPercentError12 (3)  = 0.77
  sim_focalPercentError12 (4)  = 1.12
  sim_focalPercentError12 (5)  = 0.99
  sim_focalPercentError12 (6)  = 0.41
  sim_focalPercentError12 (7)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError12 (8)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError12 (9)  = 0.01    ! not actual value, temporary placeholder
  sim_focalPercentError12 (10) = 0.01    ! not actual value, temporary placeholder
!
!
!     ...The same but for the power deposition errors.
!
!
  sim_powerPercentError1  (1)  = 0.99
  sim_powerPercentError1  (2)  = 0.46
  sim_powerPercentError1  (3)  = 0.11
  sim_powerPercentError1  (4)  = 0.02
  sim_powerPercentError1  (5)  = 0.01
  sim_powerPercentError1  (6)  = 0.01
  sim_powerPercentError1  (7)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError1  (8)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError1  (9)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError1  (10) = 0.01    ! not actual value, temporary placeholder

  sim_powerPercentError12 (1)  = 1.24
  sim_powerPercentError12 (2)  = 0.04
  sim_powerPercentError12 (3)  = 0.10
  sim_powerPercentError12 (4)  = 0.05
  sim_powerPercentError12 (5)  = 0.09
  sim_powerPercentError12 (6)  = 0.03
  sim_powerPercentError12 (7)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError12 (8)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError12 (9)  = 0.01    ! not actual value, temporary placeholder
  sim_powerPercentError12 (10) = 0.01    ! not actual value, temporary placeholder
!
!
!     ...Convert the info from refinement level -> individual rays.
!
!
  sim_rayFexitPercentError1  = sim_focalPercentError1  (sim_refinementLevel)
  sim_rayFexitPercentError12 = sim_focalPercentError12 (sim_refinementLevel)
  sim_rayPexitPercentError1  = sim_powerPercentError1  (sim_refinementLevel)
  sim_rayPexitPercentError12 = sim_powerPercentError12 (sim_refinementLevel)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_setErrorBars
