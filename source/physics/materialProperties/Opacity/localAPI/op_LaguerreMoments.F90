!!****if* source/physics/materialProperties/Opacity/localAPI/op_LaguerreMoments
!!
!! NAME
!!
!!  op_LaguerreMoments
!!
!! SYNOPSIS
!!
!!  call op_LaguerreMoments (integer (in)    :: n,
!!                           integer (in)    :: nLagT,
!!                           real    (in)    :: alpha,
!!                           real    (in)    :: beta,
!!                           real    (in)    :: T,
!!                           real    (inout) :: LagT (1:nLagT),
!!                           real    (out)   :: MomZero,
!!                           real    (out)   :: Mom  (1:n))
!!
!! DESCRIPTION
!!
!!  Generates up to n-th order (normalized) moments over generalized Laguerre polynomials
!!  Lag (alpha,x) and Laguerre weights x^beta * exp (-x) over the integration range {x = 0 to T}:
!!
!!                                        T
!!                                      /
!!                                     |   alpha     beta  -x
!!             Mom (T,alpha,beta) =    |  Lag   (x) x     e   dx
!!                i                    |     i
!!                                    /
!!                                   0
!!
!!  by using the following recursion:
!!
!!                                                     alpha     -T  beta+1
!!               Mom    = (beta - alpha - i) Mom   -  Lag   (T) e   T            i = 0,1,...,n-1
!!                  i+1                         i        i
!!
!!
!!  with starting value Mom  = MomZero.
!!                         0
!!
!!  The obtained moments are normalized by dividing all the moments by 1/MomZero . The 0-th
!!  moment is evaluated in this routine.
!!
!!  Stable form of moment recursion
!!  -------------------------------
!!
!!  The recursion is very stable in the forward direction for a wide variety of alpha, beta and
!!  T ranges. Special care has to be taken however to stay within computational upper and lower
!!  numerical bounds for large i indices. The following 'tricks' are applied:
!!
!!
!!      1) Use rescaled instead of standard generalized Laguerre polynomials (see comments in
!!         the routine evaluating these polynomials.
!!
!!      2) Move all of the factors into exponential form and perform only one exponential
!!         call. This avoids numerical problems with evaluation of products between very
!!         small and/or very large numbers in which the individual terms might run out of
!!         computational bounds but the product does not (for example, e^-T T^beta for large
!!         values of T).
!!
!!
!!  Point 1) requires introduction of the Laguerre polynomial rescaling factor into the
!!  above moment recursion form. The technique of point 2) is shown in what follows.
!!
!!  The introduction of the Laguerre polynomial rescaling factor:
!!
!!                             k=i
!!                 P (T)  =   Prod  (T +/- 2k)
!!                  i          k=1
!!
!!
!!  means that:
!!
!!               alpha           1   alpha
!!              Lag     (T)  =  --- Lag     (T)
!!                 i,res         P     i
!!                                i
!!
!!  Our moment recursion becomes:
!!
!!                                           alpha             -T  beta+1
!!     Mom    = (beta - alpha - i) Mom   -  Lag     (T) P (T) e   T
!!        i+1                         i        i,res     i
!!
!!
!!  The last term can now be evaluated using the strategy mentioned in point 2) above.
!!  We have:
!!
!!
!!   alpha             -T  beta+1       alpha                            -T  beta+1
!!  Lag     (T) P (T) e   T        =   Lag     (T) sign (P (T)) |P (T)| e   T
!!     i,res     i                        i,res           i       i
!!
!!
!!                                alpha                     (-T + [beta+1] * ln T + ln |P (T)|)
!!                           =   Lag     (T) sign (P (T)) e                              i
!!                                  i,res           i
!!
!!                                                                                   i
!!                                alpha                     (-T + [beta+1] * ln T + Sum ln |T+/-2k|)
!!                           =   Lag     (T) sign (P (T)) e                         k=1
!!                                  i,res           i
!!
!!                                                                                             i
!!                                                         (-T + [beta+1] * ln T + ln |Lag| + Sum ln |T+/-2k|)
!!                           =   sign (Lag) sign (P (T)) e                                    k=1
!!                                                 i
!!
!!
!!  The accumulation of the sum of the logarithmic rescaling terms as well as determination of
!!  the rescaling factor sign is done as an update at each recursion step. Note, that we trade
!!  computation time for stability. The penalty we pay for having a stable recursion is of course
!!  the repeated calls to the logarithm intrinsic function.
!!
!!
!! ARGUMENTS
!!
!!  n       : the maximum order of the moments wanted
!!  nLagT   : the declared dimension for the array holding the Laguerre polynomials
!!  alpha   : the generalizing Laguerre polynomial factor
!!  beta    : the exponential for the generalized Laguerre weight function
!!  T       : the upper integration limit for the moments (must be > 0)
!!  LagT    : will hold the T-valued rescaled generalized Laguerre polynomials for i = 1,...,n-1
!!  MomZero : the 0-th moment
!!  Mom     : the normalized moments for i = 1,...,n
!!
!!***
subroutine op_LaguerreMoments (n, nLagT, alpha, beta, T, LagT, MomZero, Mom)

  implicit none

  integer, intent (in)    :: n
  integer, intent (in)    :: nLagT
  real,    intent (in)    :: alpha
  real,    intent (in)    :: beta
  real,    intent (in)    :: T
  real,    intent (inout) :: LagT (1:nLagT)
  real,    intent (out)   :: MomZero
  real,    intent (out)   :: Mom  (1:n)

  LagT    = 0.0
  Mom     = 0.0
  MomZero = 0.0

  return
end subroutine op_LaguerreMoments
