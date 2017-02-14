!!****if* source/Simulation/SimulationMain/Cellular/sim_ranmar
!!
!! NAME
!!
!!  sim_ranmar
!!
!!
!! SYNOPSIS
!!
!!  sim_ranmar(integer(INOUT) :: ijkl,
!!               real(INOUT)  :: rvec,
!!               integer(IN)  :: len )
!!
!! DESCRIPTION
!!
!!  marsaglia and zaman random number generator. 
!!  period is 2**43 with 900 million different sequences
!!  the state of the generator (for restarts) is in the common block. 
!!
!! ARGUMENTS
!!
!!  ijkl - seed.  Set to a negative number on input to initialize or reinitialize.
!!         On output, it is set to the absolute value of input
!!  rvec - array of random numbers between 0.0 and 1.0
!!  len  - length of rvec
!!
!! PARAMETERS
!!
!!
!! NOTES
!!     ijkl = 54217137 is the same as the four seeds 
!!            i=12, j=34, k=56, l=78 of marsagalia and zaman 
!!
!!
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!  
!!  
!!***

subroutine sim_ranmar(ijkl,rvec,len) 
  implicit none 

  !..declarations 

  !     Parameter
  integer, INTENT(in)   ::  len 
  integer,INTENT(inout) ::  ijkl
  real,INTENT(inout)    ::  rvec(len)

  !     Local
  integer         ::  ivec,ij,kl,i,j,k,l,ii,jj,m,iff  
  integer, save   :: i97, j97
  real            ::       uni,s,t  
  real, save      :: u(97), c, cd, cm


  !  common /ranset1/ u,c,cd,cm,i97,j97 

  !..initalization initialize even if the user forgot to set ijkl to
  !..a negative value the first time.
  data             iff/0/ 

  !------------------------------------------------------------------------------

  if (ijkl .lt. 0 .or. iff .eq. 0) then 
     iff  = 1 

     if (ijkl .lt. 0) ijkl = -ijkl 

     ij   = ijkl/30082 
     kl   = ijkl - 30082 * ij 
     i    = mod(ij/177,177) + 2 
     j    = mod(ij,177) + 2 
     k    = mod(kl/169,178) + 1 
     l    = mod(kl,169) 

     do ii=1,97 
        s = 0.0
        t = 0.5
        do jj =1,24 
           m = mod(mod(i*j,179)*k,179) 
           i = j 
           j = k 
           k = m 
           l = mod(53*l + 1, 169) 
           if (mod(l*m,64) .ge. 32) s = s + t 
           t = 0.5 * t 
        enddo
        u(ii) = s 
     enddo
     c   = 362436.0/16777216.0
     cd  = 7654321.0/16777216.0
     cm  = 16777213.0/16777216.0
     i97 = 97 
     j97 = 33 
  end if

  !..normal execution starts here
  do ivec=1,len 
     uni = u(i97) - u(j97) 
     if (uni .lt. 0.0) uni = uni + 1.0
     u(i97) = uni 
     i97 = i97 - 1 
     if (i97 .eq. 0) i97 = 97 
     j97 = j97 - 1 
     if (j97 .eq. 0) j97 = 97 
     c = c - cd 
     if (c .lt. 0.0) c = c + cm 
     uni = uni - c 
     if (uni .lt. 0.0) uni = uni + 1.0
     rvec(ivec) = uni 
  enddo

  return 

end subroutine sim_ranmar
