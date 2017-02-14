Module fl_fsConeData
  implicit none

  integer, save :: nc12it, nne22it, nldent
  integer, save :: ilt, iut, jlt, jut, klt, kut
  integer, save :: lda, lna, lca !! cached indices for ut_hunt

! convention i --> log(rho), j --> X22, k --> X12

  real, save, dimension(:), allocatable :: c12ib, ne22ib, ldenb
  real, save, dimension(:,:,:), allocatable :: stab, dstab

end Module fl_fsConeData
