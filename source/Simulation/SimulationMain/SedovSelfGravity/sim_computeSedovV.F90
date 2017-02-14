!!****if* source/Simulation/SimulationMain/SedovSelfGravity/sim_computeSedovV
!!
!! NAME
!!
!!  sim_computeSedovV
!!
!!
!! SYNOPSIS
!!   sim_computeSedovV(real(in) :: xi,
!!                     real(in) :: gamma, 
!!                     real(in) :: nu1,
!!                     real(in) :: nu2,
!!                     real(out) :: V)
!!
!! DESCRIPTION
!!   Compute the dimensionless velocity function V(xi) in the
!!   Sedov problem using bisection.  See Landau & Lifshitz for
!!   notation.
!!
!! ARGUMENTS
!!   xi - coordinate of the point to calculate velocity
!!   gamma -  gas gamma
!!   nu1 -
!!   nu2 -
!!   V - calculated velocity
!! 
!!***

subroutine sim_computeSedovV (xi, gamma, nu1, nu2, V)

  !========================================================================
  implicit none
  real, intent(IN)::     xi, gamma, nu1, nu2
  real, intent(OUT) :: V

  real      VL, VR, xiL, xiR, Vmid, ximid, sedov_vfunc, tolerance, & 
       &            logxi
  integer   n_iter, n_iter_max
  parameter (n_iter_max = 500, tolerance = 1.E-6)
  
  !===============================================================================
  
  if (xi .le. 1.E-6) then         ! Use asymptotic xi->0 solution.
     
     V = 1./gamma
     
  else                            ! Do bisection root search.
     
     logxi = alog(xi)
     VL = 1./gamma
     VR = 2./gamma
     xiL = sedov_vfunc(VL, gamma, nu1, nu2)
     xiR = sedov_vfunc(VR, gamma, nu1, nu2)
     n_iter = 1
10   Vmid = 0.5 * (VL + VR)
     ximid = sedov_vfunc(Vmid, gamma, nu1, nu2)
     if ((abs(ximid - logxi) .le. tolerance) .or. & 
          &          (n_iter .gt. n_iter_max)) goto 20
     n_iter = n_iter + 1
     if (ximid .gt. logxi) then
        VR = Vmid
     else
        VL = Vmid
     endif
     goto 10
20   if (n_iter .gt. n_iter_max) & 
          &      write (*,*) 'sim_computeSedovV:  did not reach ', & 
          &                  'max precision for xi = ', xi
     V = Vmid
     
  endif
  
  !==================================================================

  return
end subroutine sim_computeSedovV


!**********************************************************************

!  Function:    sedov_vfunc()

!  Description: Function to use in bisection root search (sim_computeSedovV()).

real function sedov_vfunc (V, gamma, nu1, nu2)
  
  implicit none
  real,intent(IN):: V, gamma, nu1, nu2
  real gamp1, gamm1, gam7, k, xi
  
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam7  = 7. - gamma
  k     = gamp1 / gamm1
  
  xi = nu1*alog(5.-(3.*gamma-1.)*V) + & 
       nu2*alog(gamma*V-1.) - & 
       nu1*alog(gam7/gamp1) - nu2*alog(gamm1/gamp1) - & 
       2.*alog(0.5*gamp1)
  sedov_vfunc = 0.2 * xi
  return
end function sedov_vfunc

!*******************************************************************************

!  Routine:     find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

        subroutine find (x, N, x0, i)

        integer N, i, il, ir, im
        real    x(N), x0

        if (x0 .lt. x(1)) then

          i = 0

        elseif (x0 .gt. x(N)) then

          i = N

        else

          il = 1
          ir = N
10          if (ir .eq. il+1) goto 20
            im = (il + ir) / 2
            if (x(im) .gt. x0) then
              ir = im
            else
              il = im
            endif
            goto 10
20        i = il

        endif

        return
        end
