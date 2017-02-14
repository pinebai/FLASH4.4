!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/st_ounoiseupdate
!!
!! NAME
!!
!!  st_ounoiseupdate
!!
!! SYNOPSIS
!!
!!  st_ounoiseupdate(integer, intent (IN)  :: vectorlength,
!!                   real, intent (INOUT)  :: vector,
!!                   real, intent (IN)  :: variance,
!!                   real, intent (IN)  :: dt,
!!                   real, intent (IN)  :: ts)
!!
!! DESCRIPTION
!!
!! Subroutine updates a vector of real values according to an algorithm
!!   that generates an Ornstein-Uhlenbeck, or "colored noise" sequence.
!!
!! The sequence x_n is a Markov process that takes the previous value, 
!!   weights by an exponential damping factor with a given correlation
!!   time "ts", and drives by adding a Gaussian random variable with
!!   variance "variance", weighted by a second damping factor, also
!!   with correlation time "ts". For a timestep of dt, this sequence 
!!   can be written as :
!!
!!     x_n+1 = f x_n + sigma * sqrt (1 - f**2.) z_n
!!
!! where f = exp (-dt / ts), z_n is a Gaussian random variable drawn
!! from a Gaussian distribution with unit variance, and sigma is the
!! desired variance of the OU sequence. (See Bartosch, 2001).
!!
!! The resulting sequence should satisfy the properties of zero mean,
!!   and stationary (independent of portion of sequence) RMS equal to
!!   "variance". Its power spectrum in the time domain can vary from  
!!   white noise to "brown" noise (P (f) = const. to 1 / f^2).
!! 
!! References :
!!    Bartosch, 2001
!! http://octopus.th.physik.uni-frankfurt.de/~bartosch/publications/IntJMP01.pdf
!!   Finch, 2004
!! http://pauillac.inria.fr/algo/csolve/ou.pdf
!!         Uhlenbeck & Ornstein, 1930
!! http://prola.aps.org/abstract/PR/v36/i5/p823_1
!!
!!
!! ARGUMENTS
!!
!!   vectorlength : length of vector to be updated
!!
!!   vector :       vector to be updated
!!
!!   variance :     variance of the distributio
!! 
!!   dt :           timestep
!!
!!   ts :           correlation time
!!
!!
!!
!!***


subroutine st_ounoiseupdate (vectorlength, vector, variance, dt, ts)
  
  implicit none
  
  real, intent (IN)                         :: variance, dt, ts
  integer, intent (IN)              :: vectorlength 
  real, intent (INOUT)              :: vector (vectorlength)
  real                              :: grnval, damping_factor
  integer                           :: i
  
  damping_factor = exp (-dt / ts)
  
  do i = 1, vectorlength
     call st_grn (grnval)
     vector (i) = vector (i) * damping_factor + variance *   &
          sqrt (1. - damping_factor**2.) * grnval
  end do
  
  return 
  
end subroutine st_ounoiseupdate

