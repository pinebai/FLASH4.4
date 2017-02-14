Module pt_dpdData
  integer, save :: pt_dpdUpdateCycle=1;
  integer, save :: pt_dpdNstep
  real, save :: pt_dpdLambda
  real, allocatable, dimension(:,:), save :: particles2
end Module pt_dpdData
