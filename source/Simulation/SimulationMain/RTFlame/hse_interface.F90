! Dean Townsley 2008

Module hse_interface

interface
subroutine sim_hse_step(dens, temp, ye, sumy, n, inputg, delta, direction, order, mode)
  implicit none
  real, dimension(:), intent(IN)    :: ye, sumy
  real, dimension(:), intent(INOUT) :: dens, temp
  integer, intent(IN)               :: n               ! index to update
  real, intent(IN)                  :: inputg, delta
  integer, intent(IN)               :: direction, order, mode
end subroutine sim_hse_step
end interface

interface
subroutine flame_hse(flam, dens, temp, ye, sumy, &
                     dens_u, temp_u, ye_u, sumy_u, &
                     dens_b, temp_b, ye_b, sumy_b, &
                     fposition, x1, dx, grav, nfill)
  implicit none
  real, dimension(:), intent(OUT)  :: flam, dens, temp, ye, sumy
  real, intent(IN)                 :: dens_u, temp_u, ye_u, sumy_u
  real, intent(IN)                 :: dens_b, temp_b, ye_b, sumy_b
  real, intent(IN)                 :: fposition, x1, dx, grav
  integer, intent(IN)              :: nfill
  ! x0 is coordinate of first element in return arrays, dx is spacing and nfill
  ! is the number of cells that this function should fill
end subroutine flame_hse
end interface

end Module hse_interface
