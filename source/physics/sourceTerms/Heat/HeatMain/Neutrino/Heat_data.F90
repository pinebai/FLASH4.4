module Heat_data
  
  real, save :: ht_Tneut, ht_Lneut
  logical, save :: useHeat
  integer, save :: ht_meshMe, ht_numProcs
  real, save :: ht_bounceTime
  logical, save :: ht_postBounce, ht_useEntr
  real, save :: ht_heatTimeFac

end module Heat_data
