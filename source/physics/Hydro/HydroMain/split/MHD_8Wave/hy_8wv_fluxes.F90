!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_fluxes
!!
!! NAME
!!
!!  hy_8wv_fluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_fluxes(real(IN)    :: Um(NUNK_VARS,n),
!!                real(IN)    :: Up(NUNK_VARs,n),
!!                real(OUT)   :: Flux(NFLUXES,n),
!!                real(OUT)   :: speed(n),
!!                real(OUT)   :: vint(n),
!!                integer(IN) :: n,
!!                integer(IN) :: dir )
!!
!!
!! DESCRIPTION
!!
!!  Given interpolated data this function computes interface MHD fluxes
!!  using Roe-type linearization of the full system.
!!  Reference: Powell et al, J. Comput. Phys., 154(2), 284, 1999.
!!
!!
!! ARGUMENTS
!!
!!  Um,Up-       Arrays of left and right interpolated variables
!!  Flux -       Array of MHD fluxes
!!  speed-       Array of fastest characteristic velocities passed
!!               to the calling routine that computes CFL condition
!!  vint -       Array of interface velocities passed to the calling
!!               routine that advances species in Lagrangian fashion
!!  n    -       Size of arrays in the sweep direction
!!  dir  -       Sweep direction
!!
!!***

subroutine hy_8wv_fluxes(Um,Up,Flux,speed,vint,n,dir)

  implicit none

#include "Flash.h"
#include "constants.h"

  !!$ Argument list -------------------------------------
  integer, INTENT(in) :: n
  real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: Um,Up
  real, DIMENSION(NFLUXES,n), INTENT(out) :: Flux
  real, DIMENSION(n), INTENT(out) :: speed, vint
  integer, INTENT(in):: dir
  !!$ ---------------------------------------------------

  real, PARAMETER :: sqhalf = 0.7071067811865475244

  real :: s,d,c,c2,un,cs,cf,as,af,bn,sbn
  real :: csl,csr,cfl,cfr,pl,pr,el,er,ed,ep
  real :: sqrtd,B2,V2
  integer :: i,VELN_VAR,MOMN_FLUX,MAGN_VAR,MAGN_FLUX

  real, DIMENSION(3) :: V,B
  real, DIMENSION(3) :: Bt,t

  logical :: RoeAvg = .true.

  Flux  = 0.
  speed = 0.
  vint  = 0.


  select case (dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     MAGN_VAR = MAGX_VAR
     
     MOMN_FLUX = XMOM_FLUX
     MAGN_FLUX = MAGX_FLUX

  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     MAGN_VAR = MAGY_VAR
     
     MOMN_FLUX = YMOM_FLUX
     MAGN_FLUX = MAGY_FLUX

  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     MAGN_VAR = MAGZ_VAR
     
     MOMN_FLUX = ZMOM_FLUX
     MAGN_FLUX = MAGZ_FLUX
  end select

  !        Up(i-1)=UL  UR=Um(i)
  !                |   |
  ! -------*---------|---------*---------
  !       i-1       i-1/2      i
  !
  do i = 3, n-1

     ! Total pressure and energy
     B2 = 0.5*(Up(MAGX_VAR,i-1)**2+Up(MAGY_VAR,i-1)**2+Up(MAGZ_VAR,i-1)**2)
     pl = Up(PRES_VAR,i-1)+B2
     el = Up(DENS_VAR,i-1)*(Up(EINT_VAR,i-1)+0.5*(Up(VELX_VAR,i-1)**2+Up(VELY_VAR,i-1)**2+Up(VELZ_VAR,i-1)**2))+B2

     B2 = 0.5*(Um(MAGX_VAR, i )**2+Um(MAGY_VAR, i )**2+Um(MAGZ_VAR, i )**2)
     pr = Um(PRES_VAR, i )+B2
     er = Um(DENS_VAR, i )*(Um(EINT_VAR, i )+0.5*(Um(VELX_VAR, i )**2+Um(VELY_VAR, i )**2+Um(VELZ_VAR, i )**2))+B2


     ! Roe-average variables
     s = sqrt(Um(DENS_VAR,i)/Up(DENS_VAR,i-1)) ! sqrt(rho_R/rho_L)
     d = s*Up(DENS_VAR,i-1)                    ! rho^hat = sqrt(rho_R)*sqrt(rho_L)
     sqrtd = sqrt(d)                           ! sqrt(rho^hat)
     s = 1./(1.+s)                             ! sqrt(rho_L)/[sqrt(rho_L) + sqrt(rho_R)]
     c = 1.-s                                  ! sqrt(rho_R)/[sqrt(rho_L) + sqrt(rho_R)]

     ed = 0.
     ep = s*(1.+Up(DENS_VAR,i-1)*Up(EINT_VAR,i-1)/Up(PRES_VAR,i-1))/Up(GAMC_VAR,i-1)+ &
          c*(1.+Um(DENS_VAR, i )*Um(EINT_VAR, i )/Um(PRES_VAR, i ))/Um(GAMC_VAR, i )

     V  = s*Up(VELX_VAR:VELZ_VAR,i-1)+c*Um(VELX_VAR:VELZ_VAR,i)
     un = V(dir)
     V2 = V(IAXIS)**2+V(JAXIS)**2+V(KAXIS)**2
     B  = s*Um(MAGX_VAR:MAGZ_VAR,i)+c*Up(MAGX_VAR:MAGZ_VAR,i-1)
     bn = B(dir)
     B2 = B(IAXIS)**2+B(JAXIS)**2+B(KAXIS)**2
     c2 = ((s*(er+pr)+c*(el+pl)-B2)/d-0.5*V2-ed)/ep

     sbn = sign(1.,bn)

     ! Begin fast and slow speeds
     cf  = 0.5*(c2+B2/d)
     cfl = pl+0.5*(Up(GAMC_VAR,i-1)-2.)*Up(PRES_VAR,i-1)
     cfr = pr+0.5*(Um(GAMC_VAR, i )-2.)*Um(PRES_VAR, i )
     cs  = abs(cf-sqrt(abs(cf*cf-c2*bn*bn/d)))
     csl = abs((cfl-sqrt(abs(cfl*cfl-Up(GAMC_VAR,i-1)*Up(PRES_VAR,i-1)*Up(MAGN_VAR,i-1)**2))))/Up(DENS_VAR,i-1)
     csr = abs((cfr-sqrt(abs(cfr*cfr-Um(GAMC_VAR, i )*Um(PRES_VAR, i )*Um(MAGN_VAR, i )**2))))/Um(DENS_VAR, i )
     cf  = 2.*cf-cs
     cfl = 2.*cfl/Up(DENS_VAR,i-1)-csl
     cfr = 2.*cfr/Um(DENS_VAR, i )-csr

     ! Renormalization coefficients
     if( cf-cs > 0. ) then
        af = min(1.,max(0.,(c2-cs)/(cf-cs)))
        as = sqrt(1.-af)
        af = sqrt(af)
     else
        as = sqhalf
        af = sqhalf
     end if

     ! Finish fast and slow speeds
     c  = sqrt(c2)
     cs = sqrt(cs)
     cf = sqrt(cf)
     csl = sqrt(csl)
     cfl = sqrt(cfl)
     csr = sqrt(csr)
     cfr = sqrt(cfr)

     select case (dir)
     case (SWEEP_X)
        s = sqrt(B(JAXIS)**2+B(KAXIS)**2)
        if( s > 0. ) then
           Bt(JAXIS) = B(JAXIS)/s
           Bt(KAXIS) = B(KAXIS)/s
        else
           Bt(JAXIS) = sqhalf
           Bt(KAXIS) = sqhalf
        end if
        Bt(IAXIS) = 0.
        t(IAXIS) = 0.
        t(JAXIS) = -Bt(KAXIS)
        t(KAXIS) =  Bt(JAXIS)

     case (SWEEP_Y)
        s = sqrt(B(IAXIS)**2+B(KAXIS)**2)
        if( s > 0. ) then
           Bt(IAXIS) = B(IAXIS)/s
           Bt(KAXIS) = B(KAXIS)/s
        else
           Bt(IAXIS) = sqhalf
           Bt(KAXIS) = sqhalf
        endif
        Bt(JAXIS) = 0.
        t(JAXIS) = 0.
        t(IAXIS) =  Bt(KAXIS)
        t(KAXIS) = -Bt(IAXIS)

     case (SWEEP_Z)
        s = sqrt(B(IAXIS)**2+B(JAXIS)**2)
        if( s > 0. ) then
           Bt(IAXIS) = B(IAXIS)/s
           Bt(JAXIS) = B(JAXIS)/s
        else
           Bt(IAXIS) = sqhalf
           Bt(JAXIS) = sqhalf
        end if
        Bt(KAXIS) = 0.
        t(KAXIS) = 0.
        t(IAXIS) = -Bt(JAXIS)
        t(JAXIS) =  Bt(IAXIS)

     end select

     ! MHD flux
     if( un >= 0. ) then

        ! Left Flux
        Flux(DENS_FLUX,i) = Up(DENS_VAR,i-1)*Up(VELN_VAR,i-1)
        Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(DENS_FLUX,i)*Up(VELX_VAR:VELZ_VAR,i-1)-Up(MAGN_VAR,i-1)*Up(MAGX_VAR:MAGZ_VAR,i-1)
        Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)+pl
        Flux(ENER_FLUX,i) = Up(VELN_VAR,i-1)*(el+pl)-Up(MAGN_VAR,i-1)* &
                            dot_product(Up(VELX_VAR:VELZ_VAR,i-1),Up(MAGX_VAR:MAGZ_VAR,i-1))
        Flux(MAGX_FLUX:MAGZ_FLUX,i) = Up(VELN_VAR,i-1)*Up(MAGX_VAR:MAGZ_VAR,i-1)-Up(MAGN_VAR,i-1)*Up(VELX_VAR:VELZ_VAR,i-1)

        ! Left fast wave
        s = fix_minus(Up(VELN_VAR,i-1)-cfl,Um(VELN_VAR,i)-cfr,un-cf)
        if( s < 0. ) then
           s = (0.5*s/c2)*(af*(Um(PRES_VAR,i)-Up(PRES_VAR,i-1)-d*cf*(Um(VELN_VAR,i)-Up(VELN_VAR,i-1)))+ &
                as*dot_product(Bt,(c*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1))+ &
                (d*cs*sbn)*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))))
           Flux(DENS_FLUX,i) = Flux(DENS_FLUX,i)+s*af
           Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)-s*af*cf
           Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)+(s*af)*V+(s*as*cs*sbn)*Bt
           Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)+(s*as*c/sqrtd)*Bt
           Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+s*(af*(0.5*V2+ed-cf*un+ep*c2)+ &
                as*dot_product(Bt,(c/sqrtd)*B+(cs*sbn)*V))

           ! Left Alfven wave
           s = un-sbn*bn/sqrtd
           if( s < 0. ) then
              s = (0.5*s)*dot_product(t,d*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))+ &
                   (sbn*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1)))
              Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)+s*t
              Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)+(s*sbn/sqrtd)*t
              Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+s*dot_product(t,V+(sbn/sqrtd)*B)

              ! Left slow wave
              s = fix_minus(Up(VELN_VAR,i-1)-csl,Um(VELN_VAR,i)-csr,un-cs)
              if( s < 0. ) then
                 s = (0.5*s/c2)*(as*(Um(PRES_VAR,i)-Up(PRES_VAR,i-1)-d*cs*(Um(VELN_VAR,i)-Up(VELN_VAR,i-1)))- &
                      af*dot_product(Bt,(c*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1))+ &
                      (d*cf*sbn)*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))))
                 Flux(DENS_FLUX,i) = Flux(DENS_FLUX,i)+s*as
                 Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)-s*as*cs
                 Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)+(s*as)*V-(s*af*cf*sbn)*Bt
                 Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)-(s*af*c/sqrtd)*Bt
                 Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)+s*(as*(0.5*V2+ed-cs*un+ep*c2)- &
                      af*dot_product(Bt,(c/sqrtd)*B+(cf*sbn)*V))

              end if
           end if
        end if

     else

        ! Right flux
        Flux(DENS_FLUX,i) = Um(DENS_VAR,i)*Um(VELN_VAR,i)
        Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(DENS_FLUX,i)*Um(VELX_VAR:VELZ_VAR,i)-Um(MAGN_VAR,i)*Um(MAGX_VAR:MAGZ_VAR,i)
        Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)+pr
        Flux(ENER_FLUX,i) = Um(VELN_VAR,i)*(er+pr)-Um(MAGN_VAR,i)*dot_product(Um(VELX_VAR:VELZ_VAR,i),Um(MAGX_VAR:MAGZ_VAR,i))
        Flux(MAGX_FLUX:MAGZ_FLUX,i) = Um(VELN_VAR,i)*Um(MAGX_VAR:MAGZ_VAR,i)-Um(MAGN_VAR,i)*Um(VELX_VAR:VELZ_VAR,i)

        ! Right fast wave
        s = fix_plus(Up(VELN_VAR,i-1)+cfl,Um(VELN_VAR,i)+cfr,un+cf)
        if( s > 0. ) then
           s = (0.5*s/c2)*(af*(Um(PRES_VAR,i)-Up(PRES_VAR,i-1)+d*cf*(Um(VELN_VAR,i)-Up(VELN_VAR,i-1)))+ &
                as*dot_product(Bt,(c*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1))- &
                (d*cs*sbn)*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))))
           Flux(DENS_FLUX,i) = Flux(DENS_FLUX,i)-s*af
           Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)-s*af*cf
           Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)-(s*af)*V+(s*as*cs*sbn)*Bt
           Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)-(s*as*c/sqrtd)*Bt
           Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-s*(af*(0.5*V2+ed+cf*un+ep*c2)+ &
                as*dot_product(Bt,(c/sqrtd)*B-(cs*sbn)*V))

           ! Right Alfven wave
           s = un+sbn*bn/sqrtd
           if( s > 0. ) then
              s = (0.5*s)*dot_product(t,d*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))- &
                   (sbn*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1)))
              Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)-s*t
              Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)+(s*sbn/sqrtd)*t
              Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-s*dot_product(t,V-(sbn/sqrtd)*B)

              ! Right slow wave
              s = fix_plus(Up(VELN_VAR,i-1)+csl,Um(VELN_VAR,i)+csr,un+cs)
              if( s > 0. ) then
                 s = (0.5*s/c2)*(as*(Um(PRES_VAR,i)-Up(PRES_VAR,i-1)+d*cs*(Um(VELN_VAR,i)-Up(VELN_VAR,i-1)))- &
                      af*dot_product(Bt,(c*sqrtd)*(Um(MAGX_VAR:MAGZ_VAR,i)-Up(MAGX_VAR:MAGZ_VAR,i-1))- &
                      (d*cf*sbn)*(Um(VELX_VAR:VELZ_VAR,i)-Up(VELX_VAR:VELZ_VAR,i-1))))
                 Flux(DENS_FLUX,i) = Flux(DENS_FLUX,i)-s*as
                 Flux(MOMN_FLUX,i) = Flux(MOMN_FLUX,i)-s*as*cs
                 Flux(XMOM_FLUX:ZMOM_FLUX,i) = Flux(XMOM_FLUX:ZMOM_FLUX,i)-(s*as)*V-(s*af*cf*sbn)*Bt
                 Flux(MAGX_FLUX:MAGZ_FLUX,i) = Flux(MAGX_FLUX:MAGZ_FLUX,i)+(s*af*c/sqrtd)*Bt
                 Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-s*(as*(0.5*V2+ed+cs*un+ep*c2)- &
                      af*dot_product(Bt,(c/sqrtd)*B-(cf*sbn)*V))

              end if
           end if
        end if
     end if

     s = 0.5*abs(un)*(Um(MAGN_VAR,i)-Up(MAGN_VAR,i-1))
     Flux(MAGN_FLUX,i) = Flux(MAGN_FLUX,i)-s
     Flux(ENER_FLUX,i) = Flux(ENER_FLUX,i)-s*bn

     speed(i) = abs(un)+cf
     vint(i)  = un

  end do

  vint(2) = Um(VELN_VAR,2)
  vint(n) = Up(VELN_VAR,n-1)

contains

  function fix_plus(al,ar,a)
    real, intent(in) :: al,ar,a
    real :: d,fix_plus

    d = max(0.,4.*(ar-al))

    if( abs(a) >= 0.5*d ) then
       fix_plus = max(0.,a)
    else
       fix_plus = 0.5*(a*a/d+a+0.25*d)
    end if
  end function fix_plus

  function fix_minus(al,ar,a)
    real, intent(in) :: al,ar,a
    real :: d,fix_minus

    d = max(0.,4.*(ar-al))

    if( abs(a) >= 0.5*d ) then
       fix_minus = min(0.,a)
    else
       fix_minus = -0.5*(a*a/d-a+0.25*d)
    end if
  end function fix_minus

end subroutine hy_8wv_fluxes
