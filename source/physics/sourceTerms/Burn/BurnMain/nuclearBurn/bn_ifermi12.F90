!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_ifermi12
!!
!! NAME
!!
!!   bn_ifermi12
!!
!! SYNOPSIS
!!
!!   real function bn_ifermi12(real(IN) :: f)
!!
!! ARGUMENTS
!!
!!  f -- umm, sorry don't know what f is in ifermi12
!!  
!! DESCRIPTION
!!
!!  routine ifermi12 does an inverse fermi integral for bn_sneutx
!!  this routine applies a rational function expansion to get the inverse
!!  fermi-dirac integral of order 1/2 when it is equal to f.
!!  maximum error is 4.19e-9.   reference: antia apjs 84,101 1993
!!
!!***

real function bn_ifermi12(f)

  use bn_interface, ONLY:  bn_sneutx

  implicit none

  !!  declare arguments
  real, intent(IN) :: f

  !! declare local variables
  integer          :: i
  integer, save    :: m1,k1,m2,k2
  real             :: rn,den,ff
  real, save       :: an,a1(12),b1(12),a2(12),b2(12)


  !!  load the coefficients of the expansion
  data  an,m1,k1,m2,k2 /0.5e0, 4, 3, 6, 5/
  data  (a1(i),i=1,5)/ 1.999266880833e4,   5.702479099336e3, & 
       &     6.610132843877e2,   3.818838129486e1, & 
       &     1.0e0/
  data  (b1(i),i=1,4)/ 1.771804140488e4,  -2.014785161019e3, & 
       &     9.130355392717e1,  -1.670718177489e0/
  data  (a2(i),i=1,7)/-1.277060388085e-2,  7.187946804945e-2,  & 
       &                    -4.262314235106e-1,  4.997559426872e-1, & 
       &                    -1.285579118012e0,  -3.930805454272e-1, & 
       &     1.0e0/
  data  (b2(i),i=1,6)/-9.745794806288e-3,  5.485432756838e-2, & 
       &                    -3.299466243260e-1,  4.077841975923e-1, & 
       &                    -1.145531476975e0,  -6.067091689181e-2/


  if (f .lt. 4.0e0) then
     rn = f + a1(m1)
     do i=m1-1,1,-1
        rn = rn*f + a1(i)
     enddo
     den = b1(k1+1)
     do i=k1,1,-1
        den = den*f + b1(i)
     enddo
     bn_ifermi12 = log(f * rn/den)

  else
     ff = 1.0e0/f**(1.0e0/(1.0e0 + an))
     rn = ff + a2(m2)
     do i=m2-1,1,-1
        rn = rn*ff + a2(i)
     enddo
     den = b2(k2+1)
     do i=k2,1,-1
        den = den*ff + b2(i)
     enddo
     bn_ifermi12 = rn/(den*ff)
  end if

  return

end function bn_ifermi12



