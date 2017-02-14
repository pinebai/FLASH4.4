!!****if* source/Simulation/SimulationMain/RTFlame/flame_hse
!!
!! NAME
!!
!!  flame_hse
!!
!! SYNOPSIS
!!
!!  call flame_hse(real, dimension(:), intent(OUT)  :: flam,
!!                 real, dimension(:), intent(OUT)  :: dens,
!!                 real, dimension(:), intent(OUT)  :: temp,
!!                 real, dimension(:), intent(OUT)  :: ye,
!!                 real, dimension(:), intent(OUT)  :: sumy,
!!                 real, intent(IN)  :: dens_u,
!!                 real, intent(IN)  :: temp_u,
!!                 real, intent(IN)  :: ye_u,
!!                 real, intent(IN)  :: sumy_u,
!!                 real, intent(IN)  :: dens_b,
!!                 real, intent(IN)  :: temp_b,
!!                 real, intent(IN)  :: ye_b,
!!                 real, intent(IN)  :: sumy_b,
!!                 real, intent(IN)  :: fposition,
!!                 real, intent(IN)  :: x1,
!!                 real, intent(IN)  :: dx,
!!                 real, intent(IN)  :: grav,
!!                 integer, intent(IN)  :: nfill)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   flam : 
!!
!!   dens : 
!!
!!   temp : 
!!
!!   ye : 
!!
!!   sumy : 
!!
!!   dens_u : 
!!
!!   temp_u : 
!!
!!   ye_u : 
!!
!!   sumy_u : 
!!
!!   dens_b : 
!!
!!   temp_b : 
!!
!!   ye_b : 
!!
!!   sumy_b : 
!!
!!   fposition : 
!!
!!   x1 : 
!!
!!   dx : 
!!
!!   grav : 
!!
!!   nfill : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Dean Townsley 2008
!
! initialize flame located at fposition in HSE
subroutine flame_hse(flam, dens, temp, ye, sumy, &
                     dens_u, temp_u, ye_u, sumy_u, &
                     dens_b, temp_b, ye_b, sumy_b, &
                     fposition, x1, dx, grav, nfill)
  use hse_interface, ONLY: sim_hse_step
  use Flame_interface, ONLY : Flame_getProfile, Flame_getWidth
  use Simulation_data, ONLY: HSE_FORWARD, HSE_BACKWARD, HSE_SETTEMP
  implicit none
  real, dimension(:), intent(OUT)  :: flam, dens, temp, ye, sumy
  real, intent(IN)                 :: dens_u, temp_u, ye_u, sumy_u
  real, intent(IN)                 :: dens_b, temp_b, ye_b, sumy_b
  real, intent(IN)                 :: fposition, x1, dx, grav
  integer, intent(IN)              :: nfill
  ! x1 is (center) coordinate of first element in return arrays, dx is spacing and nfill
  ! is the number of cells that this function should fill

  real :: flamewidth, x
  integer :: i, refi


!=================================================================

  ! the flame width tells us the size of the region above and below
  ! the flame position where flam is something other than 0 or 1
  call Flame_getWidth(flamewidth)
  ! set up flame profile first without worrying about HSE
  do i = 1, nfill
     x = x1+dx*(i-1)
     if (x > fposition+flamewidth) then
        flam(i) = 0.0
     else if (x < fposition-flamewidth) then
        flam(i) = 1.0
     else
        call Flame_getProfile(x-fposition, flam(i))
     endif
     ! set up temp and material properties
     temp(i) = (1.0-flam(i))*temp_u + flam(i)*temp_b
     ye(i)   = (1.0-flam(i))*ye_u   + flam(i)*ye_b
     sumy(i) = (1.0-flam(i))*sumy_u + flam(i)*sumy_b
  enddo
     
  ! the reference point is just ahead of (above) flame edge
  refi = max(1,int( (fposition+flamewidth-x1)/dx + 1 ) )
  dens(refi) = dens_u
  !! first go up from reference point
  ! guess value
  dens(refi+1) = dens(refi)
  ! first step is first order
  call sim_hse_step(dens,temp,ye,sumy,refi+1,grav,dx,HSE_FORWARD,1,HSE_SETTEMP)
  do i = refi + 2, nfill
     dens(i) = dens(i-1)
     call sim_hse_step(dens,temp,ye,sumy,i,grav,dx,HSE_FORWARD,2,HSE_SETTEMP)
  enddo
  !! now go down from reference point
  do i = refi-1, 1, -1
     dens(i) = dens(i+1)
     call sim_hse_step(dens,temp,ye,sumy,i,grav,dx,HSE_BACKWARD,2,HSE_SETTEMP)
  enddo

  return

end subroutine
