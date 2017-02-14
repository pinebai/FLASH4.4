
! Aaron Jackson, Alan Calder

Module fl_fsAtwoodData
  implicit none

  integer, save :: fl_fsAtwoodTabLdensNum
  integer, save :: fl_fsAtwoodTabXc12Num
  real, save    :: fl_fsAtwoodTabMinLdens, fl_fsAtwoodTabMaxLdens, fl_fsAtwoodTabDLdens
  real, save    :: fl_fsAtwoodTabMinXc12, fl_fsAtwoodTabMaxXc12, fl_fsAtwoodTabDXc12
  real, save, dimension(:,:), allocatable :: fl_fsAtwoodTabA

end Module fl_fsAtwoodData
