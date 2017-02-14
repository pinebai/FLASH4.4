!!****ih* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CONe/fl_fsConeInterp
!!
!! NAME
!!
!!    fl_fsConeTable(dynamic)
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!   Manages the table for laminar flame speed look-up through X12, X22,
!!   LOG(den) from Chamulak, Brown, & Timmes.
!!  
!! ARGUMENTS
!!  
!!
!!***
!  convention i --> log(rho), j --> X22, k --> X12
subroutine fl_fsConeInterp(c12i,ne22i,den,s,ds)

  use Driver_interface, ONLY : Driver_abortFlash
  use fl_fsConeData, ONLY : ilt, iut, jlt, jut, klt, kut, &
                            nc12it, nne22it, nldent, lca, lna, lda, &
                            c12ib, ne22ib, ldenb, stab, dstab
  use ut_interpolationInterface, ONLY : ut_hunt
    
  implicit none

  real, intent(IN) :: c12i, ne22i, den
  real, intent(OUT) :: s, ds

  integer :: i, j, k

  real :: rc12, rne22, rlden
  real :: s111, s211, s121, s112, s221, s212, s122, s222
  real :: ds111, ds211, ds121, ds112, ds221, ds212, ds122, ds222

  real :: c1, c2, c3

!! find the location in the c12, ne22, and den arrays 
!! remember that we use the logs of density

  rc12  = c12i
  rne22 = ne22i
  rlden = log10(den)

  call ut_hunt(c12ib,nc12it,rc12,lca)
  call ut_hunt(ne22ib, nne22it,rne22,lna)
  call ut_hunt(ldenb,nldent,rlden,lda)

  if (lca.eq.0) then !! fail gracefully
     lca = klt
     rc12 = 0.0
  else if (lca.eq.nc12it) then
     lca = kut - 1
     rc12 = 1.0
  endif

  if (lna.eq.0) then !! fail gracefully
     lna = jlt
     rne22 = 0.0
  else if (lna.eq.nne22it) then !! warn of extrapolation
     if (ne22i > (1.0-rc12)) rne22 = 1.0-rc12
     write (6,*) '[fl_fsConeInterp] warning: Extrapolating high on ne22'
     lna = jut - 1
  endif

  if (lda.eq.0) then !! warn
     lda = ilt
!     write (6,*) '[fl_fsConeInterp] warning: Extrapolating low on log(dens)'
  else if (lda.eq.nldent) then
     lda = iut - 1
     write (6,*) '[fl_fsConeInterp] warning: Extrapolating high on log(dens)'
  endif

  if (rlden < 8.0) then !! using Timmes & Woosley data w/o Ne22
!     write (6,*) '[fl_fsConeInterp] warning: Interpolating using data w/o Ne22'
  endif

  if (lca > kut-(jut-jlt+1)) then
     write (6,*) '[fl_fsConeInterp] warning: using high-end of c12'
  endif

!! calculate the coefficients

  c3 = ( rc12 - c12ib(lca) )/              &
        ( c12ib(lca+1) - c12ib(lca) )

  c2 = ( rne22 - ne22ib(lna) )/             &
        ( ne22ib(lna+1) - ne22ib(lna) )

  c1 = ( rlden - ldenb(lda) )/             &
        ( ldenb(lda+1) - ldenb(lda))

!! build the local cubes

   s111 =  stab(lda,lna,lca)
  ds111 = dstab(lda,lna,lca)
   s211 =  stab(lda+1,lna,lca)
  ds211 = dstab(lda+1,lna,lca)
   s121 =  stab(lda,lna+1,lca)
  ds121 = dstab(lda,lna+1,lca)
   s112 =  stab(lda,lna,lca+1)
  ds112 = dstab(lda,lna,lca+1)
   s221 =  stab(lda+1,lna+1,lca)
  ds221 = dstab(lda+1,lna+1,lca)
   s212 =  stab(lda+1,lna,lca+1)
  ds212 = dstab(lda+1,lna,lca+1)
   s122 =  stab(lda,lna+1,lca+1)
  ds122 = dstab(lda,lna+1,lca+1)
   s222 =  stab(lda+1,lna+1,lca+1)
  ds222 = dstab(lda+1,lna+1,lca+1)

!! now interpolate

  s  =  exp(  ( c3                            &
                      *( (1.0-c1)*(1.0-c2)    &
                                     *s112    &
                        +      c1*(1.0-c2)    &
                                     *s212    &
                        +      c2*(1.0-c1)    &
                                     *s122    &
                        +      c1*c2          &
                                     *s222    &
                       )                      &
                +(1.0-c3)                     &
                      *( (1.0-c1)*(1.0-c2)    &
                                     *s111    &
                        +      c1*(1.0-c2)    &
                                     *s211    &
                        +      c2*(1.0-c1)    &
                                     *s121    &
                        +      c1*c2          &
                                     *s221    &
                       )                      &
              ) *log(10.0e0)                  &
           ) 

  ds  =  exp(  ( c3                           &
                      *( (1.0-c1)*(1.0-c2)    &
                                    *ds112    &
                        +      c1*(1.0-c2)    &
                                    *ds212    &
                        +      c2*(1.0-c1)    &
                                    *ds122    &
                        +      c1*c2          &
                                    *ds222    &
                       )                      &
                +(1.0-c3)                     &
                      *( (1.0-c1)*(1.0-c2)    &
                                    *ds111    &
                        +      c1*(1.0-c2)    &
                                    *ds211    &
                        +      c2*(1.0-c1)    &
                                    *ds121    &
                        +      c1*c2          &
                                    *ds221    &
                       )                      &
              ) *log(10.0e0)                  &
           ) 

  return
end subroutine
