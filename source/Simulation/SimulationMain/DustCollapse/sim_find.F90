!!****if* source/Simulation/SimulationMain/DustCollapse/sim_find
!!
!! NAME
!!  sim_find
!!
!! SYNOPSIS
!!
!!  sim_find(real(N), intent(IN) :: x,
!!           integer, intent(IN) :: N,
!!           real,    intent(IN) :: x0,
!!           integer, intent(OUT):: i)
!!
!! DESCRIPTION
!!
!!      Given a monotonically increasing table x(N) and a test value
!!      x0, return the index i of the largest table value less than
!!      or equal to x0 (or 0 if x0 < x(1)).  Use binary search.
!!
!! ARGUMENTS
!!   x - the table to be searched
!!   N - size of the table
!!   x0 - the test value
!!   i  - index of the largest table value < x0
!!
!!***

subroutine sim_find (x, N, x0, i)
  implicit none
  integer, intent(IN) :: N
  real,    intent(IN) :: x0
  real,dimension(N),intent(IN) :: x
  integer, intent(OUT) :: i

  integer il, ir, im
  
  if (x0 .lt. x(1)) then
     
     i = 0
     
  elseif (x0 .gt. x(N)) then
     
     i = N
     
  else
     
     il = 1
     ir = N
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il
     
  endif
  
  return
end subroutine sim_find
