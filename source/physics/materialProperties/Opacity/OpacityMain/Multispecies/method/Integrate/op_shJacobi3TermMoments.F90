!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_shJacobi3TermMoments
!!
!! NAME
!!
!!  op_shJacobi3TermMoments
!!
!! SYNOPSIS
!!
!!  call op_shJacobi3TermMoments (integer (in)    :: n,
!!                                integer (in)    :: maxM,
!!                                real    (in)    :: p,
!!                                real    (in)    :: beta,
!!                                real    (in)    :: T,
!!                                real    (inout) :: X    (1:maxM),
!!                                real    (inout) :: Y    (1:maxM),
!!                                real    (inout) :: Sinv (1:maxM),
!!                                real    (out)   :: MomZero,
!!                                real    (out)   :: Mom  (1:n))
!!
!! DESCRIPTION
!!
!!  Generates up to n-th order normalized moments over shifted Jacobi (monic) polynomials
!!  G (p,beta+1,x) and Laguerre weights x^beta * exp (-x) over the integration range {x = 0 to 1}:
!!
!!
!!                                  1                                       1
!!                                /                                        /
!!                               |                   beta  -Tx            |   beta  -Tx
!!             Mom (p,beta) =    |  G (p,beta+1,x)  x     e    dx    /    |  x     e    dx
!!                i              |   i                                    |
!!                              /                                        /
!!                             0                                        0
!!
!!  by using the following 3-term inhomogeneous recursion for the normalized moments:
!!
!!
!!                   Mom    = R  Mom   +  S  Mom     +  U           i = 0,1,...,n-1
!!                      i+1    i    i      i    i-1      i
!!
!!
!!  where (b = short for beta):
!!
!!                 i+b+1        (i+b+1)(p-1)
!!          R   =  -----  -   ----------------              i = 0,1,2,...,n-1
!!           i       T        (2i+p-1)(2i+p+1)
!!
!!
!!                   i(i+b)(i+b+1)(i+p-1)
!!          S   =  ------------------------                 i = 0,1,2,...,n-1
!!           i     (2i+p)(2i+p-2)(2i+p-1)^2
!!
!!
!!                    e^(-T)     k=i    p-b+k-2
!!          U   =  - --------  Product ---------            i = 0,1,2,...,n-1
!!           i       T * Mom     k=1    p+i+k-1
!!                          0
!!
!!  The U value for i=0 is obtained by setting the product function equal to 1 in this
!!  case.
!!
!!
!!  Stable form of moment recursion
!!  -------------------------------
!!
!!  The recursion is very unstable in the forward direction and cannot be used in that way.
!!  After several tries, it was found that the Olver technique coupled with the LU-decomposition
!!  method works best. The main literature reference is:
!!
!!            J.Wimp, Computation with Recurrence Relations, Boston MA, 1984
!!
!!  with all the citations mentioned therein. The Olver technique is a modification of
!!  Miller's algorithm that finds minimal solutions of 3-term inhomogeneous recursions by
!!  applying the recursion backwards. In our case the recursion is reformulated as:
!!
!!
!!                   Mom    =   1/S    Mom     - (R   /S   ) Mom     - (U   /S   )
!!                      i          i+1    i+2      i+1  i+1     i+1      i+1  i+1
!!
!!
!!  and the LU-decomposition algorithm then leads to the following working equations:
!!
!!
!!           x(1) = - U(1)/S(1) - 1
!!           y(1) = - 1 / (R(1)/S(1))
!!           y(i) = 1 / [- (R(i)/S(i)) + y(i-1)/S(i-1)]      i = 2,3,...,M
!!           x(i) = - U(i)/S(i) + y(i)*x(i-1)                i = 2,3,...,M
!!         Mom(M) = - x(M)*y(M)
!!         Mom(i) = - [Mom(i+1)/S(i) + x(i)] * y(i)          i = M-1,...,K,...,1
!!
!!
!!  In these equations K denotes the number of moments we are after and M denotes the
!!  total number of moments we need to calculate to achieve mantissa accuracy for the
!!  lowest K moments when applying the downward recursion.
!!
!!  Determination of M
!!  ------------------
!!
!!  From the above mentioned book by Wimp, we have an expression for the moments in
!!  terms of quantities eta as follows:
!!
!!                   inf
!!         Mom(k) =  Sum  eta(k,i)
!!                  i=k-1
!!
!!  where
!!
!!                                       n=i
!!         eta(k,i) = - x(i+1) * y(k) * Prod y(n+1)*(1/S(n))
!!                                       n=k
!!
!!
!!  For each moment Mom(k) the above sum converges to its true value with individual
!!  magnitude of the eta's decreasing as i -> infinity. To achieve a certain accuracy
!!  for the moments it is thus only necessary to find the point at which one can
!!  safely terminate the infinite summation. Once that point is found we set M = i+1
!!  and start evaluating the downward recursion scheme for the moments. Ideally one
!!  would need to track accuracy of every moment in the range k = 1,...,K, but it is
!!  postulated that once Mom(K) is accurate so will be all lower moments. Thus we have
!!  the additional working equations:
!!
!!
!!           etaProd = - y(K)
!!        etaK (K-1) = etaProd * x(K)
!!           etaProd = etaProd * (y(i+1)/S(i))          i = K,K+1,.....,accuracy test ok
!!        etaK   (i) = etaProd * x(i)                   i = K,K+1,.....,accuracy test ok
!!
!!
!!  Which variables need to be stored in arrays? The only ones we need are those which
!!  will be required for the downward evaluation of the moments: Sinv (i) = 1 / S(i),
!!  x(i) and y(i) for i=1,...,M. Since M is not known in advance we have to test for
!!  sufficient memory available for these three arrays once M has been found. All other
!!  variables like the U(i) and the etaK(i) are accumulated and/or calculated on the
!!  fly.
!!
!!
!! ARGUMENTS
!!
!!  n       : the maximum order of the moments wanted
!!  maxM    : the highest moment order that will be tried for reaching mantissa
!!            accuracy of the n moments. This is also equal to the declared memory
!!            for the auxilliary arrays X,Y,Sinv
!!  p       : the p-parameter for the shifted Jacobi polynomials G(p,beta+1,x)
!!  beta    : the exponential for the generalized Laguerre weight function
!!  T       : the upper integration limit for the moments (must be > 0)
!!  X       : will contain the x(i) values for i = 1,2,...,M
!!  Y       : will contain the y(i) values for i = 1,2,...,M
!!  Sinv    : will contain the 1/S(i) values for i = 1,2,...,M
!!  MomZero : the 0-th moment
!!  Mom     : the normalized moments for i = 1,...,n
!!
!!***
subroutine op_shJacobi3TermMoments (n, maxM, p, beta, T, X, Y, Sinv, MomZero, Mom)

  use Driver_interface,  ONLY : Driver_abortFlash
  use op_interface,      ONLY : op_LaguerreZeroMoment
  use op_numericsData,   ONLY : zero,one,two

  implicit none

  integer, intent (in)    :: n
  integer, intent (in)    :: maxM
  real,    intent (in)    :: p
  real,    intent (in)    :: beta
  real,    intent (in)    :: T
  real,    intent (out)   :: MomZero

  real,    intent (inout) :: X    (1:maxM)
  real,    intent (inout) :: Y    (1:maxM)
  real,    intent (inout) :: Sinv (1:maxM)
  real,    intent (out)   :: Mom  (1:n)

  integer :: i,K,M

  real    :: b,R,U
  real    :: etaK,etaKsum,etaKsumPrev
  real    :: MomI,MomIp1
  real    :: Tinv
  real    :: w,ww
!
!
!   ...Check proper numerical values.
!
!
  if (beta < zero) then
      call Driver_abortFlash ('[op_shJacobi3TermMoments] ERROR: Beta is < 0')
  end if

  if (T <= zero) then
      call Driver_abortFlash ('[op_shJacobi3TermMoments] ERROR: Integration limit T =< 0')
  end if

  if (maxM <= n + n) then
      call Driver_abortFlash ('[op_shJacobi3TermMoments] ERROR: maxM < 2n')
  end if
!
!
!   ...Immediate return if no moments wanted.
!
!
  call op_LaguerreZeroMoment (beta,one,T,MomZero)

  if (n == 0) then
      return
  end if
!
!
!   ...Accumulate x(i), y(i) and Sinv(i) values for i = 1,...K.
!
!
  K = n
  b = beta
  Tinv = one / T

  R        =  (b+two) * (Tinv - (p-one) / ((p+one+two) * (p+one)))
  U        = - exp (-T) * Tinv * (p-b-one) / (MomZero * (p+one))
  Sinv (1) = ((p+two) * (p+one) * (p+one)) / ((b+two) * (b+one))
  Y    (1) = - one / (R * Sinv (1))
  X    (1) = - U * Sinv (1) - one

  do i = 2,K
     w        = real (i)
     ww       = w + w
     R        = (w+b+one) * (Tinv - (p-one) / ((ww+p+one) * (ww+p-one)))
     U        = U * ((w+p-one) * (w+p-b-two)) / ((ww+p-one) * (ww+p-two))
     Sinv (i) = ((ww+p) * (ww+p-two) * (ww+p-one)**2) / (w * (w+b) * (w+b+one) * (w+p-one))
     Y    (i) = one / (- R * Sinv (i) + Y (i-1) * Sinv (i-1))
     X    (i) = - U * Sinv (i) + X (i-1) * Y (i-1)
  end do
!
!
!   ...Start evaluation of the etaK and their summation.
!
!
  etaK        = - Y (K)
  etaKsum     = etaK * X (K)              ! contribution for etaK (K-1)
  etaKsumPrev = etaKsum
!
!
!   ...Continue accumulating x(i), y(i) and Sinv(i) values for i = K+1,...M until
!      reaching mantissa accuracy for K-th moment. If M exceeds the supplied array
!      boundary maxM, a message is printed, but the calculation continues with partially
!      converged results.
!
!
  M = 0

  do i = K+1,maxM

     w        = real (i)
     ww       = w + w
     R        = (w+b+one) * (Tinv - (p-one) / ((ww+p+one) * (ww+p-one)))
     U        = U * ((w+p-one) * (w+p-b-two)) / ((ww+p-one) * (ww+p-two))
     Sinv (i) = ((ww+p) * (ww+p-two) * (ww+p-one)**2) / (w * (w+b) * (w+b+one) * (w+p-one))
     Y    (i) = one / (- R * Sinv (i) + Y (i-1) * Sinv (i-1))
     X    (i) = - U * Sinv (i) + X (i-1) * Y (i-1)

     etaK = etaK * Y (i) * Sinv (i-1)
     etaKsum = etaKsum + etaK * X (i)           ! contribution for etaK (i-1)

     if ((etaKsum / etaKsumPrev) == one) then
          M = i
          exit
     else
          etaKsumPrev = etaKsum
     end if

  end do

  if (M == 0) then
      write (*,*) '[op_shJacobi3TermMoments] PROBLEM: Mantissa accuracy not reached! Will continue...'
      M = maxM
  end if
!
!
!   ...Downward recursion of the moments from i = M,...,1. Store only the first
!      K moments.
!
!
  MomIp1 = - X (M) * Y (M)

  do i = M-1,K+1,-1
     MomI = - (MomIp1 * Sinv (i) + X (i)) * Y (i)
     MomIp1 = MomI
  end do

  Mom (K) = - (MomIp1 * Sinv (K) + X (K)) * Y (K)

  do i = K-1,1,-1
     Mom (i) = - (Mom (i+1) * Sinv (i) + X (i)) * Y (i)
  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_shJacobi3TermMoments
