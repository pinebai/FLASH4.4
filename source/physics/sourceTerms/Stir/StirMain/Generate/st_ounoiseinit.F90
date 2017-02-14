!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/st_ounoiseinit
!!
!! NAME
!!
!!  st_ounoiseinit
!!
!! SYNOPSIS
!!
!!  st_ounoiseinit(integer,intent (IN)  :: vectorlength,
!!                 integer,intent (IN)  :: iseed,
!!                 real,intent (IN)  :: variance,
!!                 real,intent (INOUT)  :: vector)
!!
!! DESCRIPTION
!!
!!
!! Subroutine initializes a vector of real numbers to be used as a 
!!   starting point for the Ornstein-Uhlenbeck, or "colored noise"
!!   generation sequence. Note that this should be invoked once at
!!   the very start of the program; "vector" values and the random
!!   seed should be checkpointed to ensure continuity across restarts.
!! 
!! The length of the vector is specified in "vectorlength", and the 
!!   random seed to be used is passed in as "seed". The sequence
!!   is initialized using Gaussian values with variance "variance".
!!
!! Please refer to st_ounoiseupdate for further details on algorithm.
!!
!! ARGUMENTS
!!
!!   vectorlength : lenght of the vector
!!
!!   iseed :  input seed for random number generator
!!
!!   variance : variance of distribution
!!
!!   vector : storage for starting vector
!!
!!
!!
!!***


subroutine st_ounoiseinit (vectorlength, iseed, variance, vector)
  
  use Stir_data, ONLY : st_randseed
  use ut_randomInterface, ONLY : ut_randomSeed
  implicit none
  
  integer,intent (IN)               :: vectorlength, iseed 
  real,intent (IN)                          :: variance 
  real,intent (INOUT)               :: vector (vectorlength)
  real                              :: grnval
  integer                           :: i
  
  !... Initialize pseudorandom sequence with random seed.
  
  st_randseed = iseed
  
  call ut_randomSeed(ut_put=st_randseed)
  
  do i = 1, vectorlength
     call st_grn (grnval)
     vector (i) = grnval * variance
  end do
  
  return
  
end subroutine st_ounoiseinit
