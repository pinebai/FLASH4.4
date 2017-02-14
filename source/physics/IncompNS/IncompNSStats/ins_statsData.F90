
module ins_statsData

  implicit none

  logical, save :: ins_statsRestart ! Restart statistics flag.
  integer, save :: ins_statsIntervalStep  ! Timesteps between stats computation.
  integer, save :: ins_statsN             ! Ensemble number.
  real,    save :: ins_statsStartTime    ! Initial time before start of stats.

end module
