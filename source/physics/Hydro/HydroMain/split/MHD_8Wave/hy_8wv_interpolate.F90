!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_interpolate
!!
!! NAME
!!
!!  hy_8wv_interpolate
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_interpolate(real(IN)  :: U(NUNK_VARS,n),
!!                     real(OUT) :: Uc(NUNK_VARS,n),
!!                     real(OUT) :: Um(NUNK_VARS,n),
!!                     real(OUT) :: Up(NUNK_VARS,n),
!!                     real(IN)  :: grav(n),
!!                     real(IN)  :: dt,
!!                     real(IN)  :: xc(n),
!!                     real(IN)  :: dx(n),
!!                     integer   :: n,
!!                     integer   :: dir)
!!
!!
!! DESCRIPTION
!!
!!  TVD interpolation of data in space and time using Hancock method
!!
!!
!! ARGUMENTS
!!
!!  U           - Array of solution data being interpolated
!!  Uc,Um,Up    - Arrays of time-interpolated data at the center
!!                and left and right interfaces of the cell
!!  grav        - Array containing gravitation acceleration
!!                component in the direction of the sweep
!!  dt          - Time step
!!  xc          - Array of cell center coordinates
!!  dx          - Array of cell sizes
!!  n           - Size of arrays in the sweep direction
!!  dir         - Sweep direction
!!
!!***

subroutine hy_8wv_interpolate(U,Uc,Um,Up,grav,dt,xc,dx,n,dir)

  use Hydro_data, ONLY : hy_eswitch, hy_meshMe

  use Eos_interface, ONLY : Eos

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"


  !!$ Argument list -------------------------------------
  integer, INTENT(IN) :: n,dir
  real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: U
  real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: Uc,Um,Up
  real, DIMENSION(n), INTENT(IN) :: grav,xc,dx
  real, INTENT(IN) :: dt
  !!$ ---------------------------------------------------

  real, PARAMETER :: sqhalf = 0.7071067811865475244
  real, DIMENSION(NUNK_VARS) :: Ux,Ut
  real, DIMENSION(NSPECIES) :: Xm,Xp
  real :: du,dul,dur,xmsum,xpsum

  real :: d,s,c,c2,cf,cs,af,as,sbn,sqd
  real, DIMENSION(MDIM) :: Bt,t

  ! Need to declare these to call eos
  real :: temp,abar,zbar
  integer :: i,ins,VELN_VAR,MAGN_VAR
  integer :: counter
  real,dimension(EOS_NUM)::eosData
  logical,dimension(EOS_VARS+1:EOS_NUM)::eosMask=.false.
  real, dimension(NSPECIES) :: eosMf
  integer :: vecLen, specIndex
  integer :: interp_eosMode = MODE_DENS_PRES
  real, PARAMETER :: epsilon = 10E-15

  select case(dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     MAGN_VAR = MAGX_VAR
  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     MAGN_VAR = MAGY_VAR
  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     MAGN_VAR = MAGZ_VAR
  end select


  do i = 2, n-1

    ! Compute characteristic variable limited gradients

    Ut = 0.
    Ux = 0.

    ! Prepare characterisitc decomposition of A(U)

    d = U(DENS_VAR,i)

    sqd = sqrt(d)
    sbn = sign(1.,U(MAGN_VAR,i))

    c2 = U(GAMC_VAR,i)*U(PRES_VAR,i)/d
    cf = 0.5*(c2+(U(MAGX_VAR,i)**2+U(MAGY_VAR,i)**2+U(MAGZ_VAR,i)**2)/d)
    cs = abs(cf-sqrt(abs(cf*cf-c2*U(MAGN_VAR,i)**2/d)))
    cf = 2.0*cf-cs

    ! Renormalization coefficients
    if( cf > cs ) then
      af = min(1.,max(0.,(c2-cs)/(cf-cs)))
      as = sqrt(1.-af)
      af = sqrt(af)
    else
      as = sqhalf
      af = sqhalf
    end if

    ! Fast, slow and sound speeds
    cf = sqrt(cf)
    cs = sqrt(cs)
    c  = sqrt(c2)

    ! Choose tangential direction
    select case(dir)
    case (SWEEP_X)

      Bt(IAXIS) = 0.
      t(IAXIS)  = 0.

      s = sqrt(U(MAGY_VAR,i)**2+U(MAGZ_VAR,i)**2)
      if( s > epsilon ) then
        Bt(JAXIS) = U(MAGY_VAR,i)/s
        Bt(KAXIS) = U(MAGZ_VAR,i)/s
      else
        Bt(JAXIS) = sqhalf
        Bt(KAXIS) = sqhalf
      end if

      t(JAXIS) = -Bt(KAXIS)
      t(KAXIS) =  Bt(JAXIS)

    case (SWEEP_Y)

      Bt(JAXIS) = 0.
      t(JAXIS)  = 0.

      s = sqrt(U(MAGX_VAR,i)**2+U(MAGZ_VAR,i)**2)
      if( s > epsilon ) then
        Bt(IAXIS) = U(MAGX_VAR,i)/s
        Bt(KAXIS) = U(MAGZ_VAR,i)/s
      else
        Bt(IAXIS) = sqhalf
        Bt(KAXIS) = sqhalf
      endif

      t(IAXIS) =  Bt(KAXIS)
      t(KAXIS) = -Bt(IAXIS)

    case (SWEEP_Z)

      Bt(KAXIS) = 0.
      t(KAXIS)  = 0.

      s = sqrt(U(MAGX_VAR,i)**2+U(MAGY_VAR,i)**2)
      if( s > epsilon ) then
        Bt(IAXIS) = U(MAGX_VAR,i)/s
        Bt(JAXIS) = U(MAGY_VAR,i)/s
      else
        Bt(IAXIS) = sqhalf
        Bt(JAXIS) = sqhalf
      end if

      t(IAXIS) = -Bt(JAXIS)
      t(JAXIS) =  Bt(IAXIS)

    end select

    ! Multiple species
    do ins = SPECIES_BEGIN, SPECIES_END
      dul = U(ins, i )-U(ins,i-1)
      dur = U(ins,i+1)-U(ins, i )

      Ux(ins) = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))
    end do

    ! Left-going fast magnetoacoustic wave
    dul = (0.5/c2)*(af*(U(PRES_VAR,i)-U(PRES_VAR,i-1)-d*cf*(U(VELN_VAR,i)-U(VELN_VAR,i-1)))+ &
                  as*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1))+ &
                                 d*cs*sbn*(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))))
    dur = (0.5/c2)*(af*(U(PRES_VAR,i+1)-U(PRES_VAR,i)-d*cf*(U(VELN_VAR,i+1)-U(VELN_VAR,i)))+ &
                  as*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i))+ &
                                 d*cs*sbn*(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(DENS_VAR)     = Ux(DENS_VAR)+du*af
    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)+du*as*cs*sbn*Bt/d
    Ux(VELN_VAR)     = Ux(VELN_VAR)-du*af*cf/d
    Ux(PRES_VAR)     = Ux(PRES_VAR)+du*af*c2
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)+du*as*c*Bt/sqd

    ! Left-going Alfven wave
    dul = 0.5*dot_product(t,(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))+ &
                  (sbn/sqd)*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1)))
    dur = 0.5*dot_product(t,(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))+ &
                  (sbn/sqd)*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i)))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)+du*t
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)+du*sbn*sqd*t

    ! Left-going slow magnetoacoustic wave
    dul = (0.5/c2)*(as*(U(PRES_VAR,i)-U(PRES_VAR,i-1)-d*cs*(U(VELN_VAR,i)-U(VELN_VAR,i-1)))- &
                  af*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1))+ &
                                 d*cf*sbn*(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))))
    dur = (0.5/c2)*(as*(U(PRES_VAR,i+1)-U(PRES_VAR,i)-d*cs*(U(VELN_VAR,i+1)-U(VELN_VAR,i)))- &
                  af*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i))+ &
                                 d*cf*sbn*(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(DENS_VAR)     = Ux(DENS_VAR)+du*as
    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)-du*af*cf*sbn*Bt/d
    Ux(VELN_VAR)     = Ux(VELN_VAR)-du*as*cs/d
    Ux(PRES_VAR)     = Ux(PRES_VAR)+du*as*c2
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)-du*af*c*Bt/sqd

    ! Entropy wave
    dul = U(DENS_VAR,i)-U(DENS_VAR,i-1)-(U(PRES_VAR,i)-U(PRES_VAR,i-1))/c2
    dur = U(DENS_VAR,i+1)-U(DENS_VAR,i)-(U(PRES_VAR,i+1)-U(PRES_VAR,i))/c2

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(DENS_VAR) = Ux(DENS_VAR)+du

    ! Divergence wave
    dul = U(MAGN_VAR,i)-U(MAGN_VAR,i-1)
    dur = U(MAGN_VAR,i+1)-U(MAGN_VAR,i)

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(MAGN_VAR) = Ux(MAGN_VAR)+du

    ! Right-going slow magnetoacoustic wave
    dul = (0.5/c2)*(as*(U(PRES_VAR,i)-U(PRES_VAR,i-1)+d*cs*(U(VELN_VAR,i)-U(VELN_VAR,i-1)))- &
                  af*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1))- &
                                 d*cf*sbn*(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))))
    dur = (0.5/c2)*(as*(U(PRES_VAR,i+1)-U(PRES_VAR,i)+d*cs*(U(VELN_VAR,i+1)-U(VELN_VAR,i)))- &
                  af*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i))- &
                                 d*cf*sbn*(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(DENS_VAR)     = Ux(DENS_VAR)+du*as
    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)+du*af*cf*sbn*Bt/d
    Ux(VELN_VAR)     = Ux(VELN_VAR)+du*as*cs/d
    Ux(PRES_VAR)     = Ux(PRES_VAR)+du*as*c2
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)-du*af*c*Bt/sqd

    ! Right-going Alfven wave
    dul = 0.5*dot_product(t,(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))- &
                  (sbn/sqd)*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1)))
    dur = 0.5*dot_product(t,(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))- &
                  (sbn/sqd)*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i)))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)+du*t
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)-du*sbn*sqd*t

    ! Right-going fast magnetoacoustic wave
    dul = (0.5/c2)*(af*(U(PRES_VAR,i)-U(PRES_VAR,i-1)+d*cf*(U(VELN_VAR,i)-U(VELN_VAR,i-1)))+ &
                  as*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i)-U(MAGX_VAR:MAGZ_VAR,i-1))- &
                                 d*cs*sbn*(U(VELX_VAR:VELZ_VAR,i)-U(VELX_VAR:VELZ_VAR,i-1))))
    dur = (0.5/c2)*(af*(U(PRES_VAR,i+1)-U(PRES_VAR,i)+d*cf*(U(VELN_VAR,i+1)-U(VELN_VAR,i)))+ &
                  as*dot_product(Bt,c*sqd*(U(MAGX_VAR:MAGZ_VAR,i+1)-U(MAGX_VAR:MAGZ_VAR,i))- &
                                 d*cs*sbn*(U(VELX_VAR:VELZ_VAR,i+1)-U(VELX_VAR:VELZ_VAR,i))))

    du = (sign(1.,dul)+sign(1.,dur))*min(abs(dul),0.25*abs(dul+dur),abs(dur))

    Ux(DENS_VAR)     = Ux(DENS_VAR)+du*af
    Ux(VELX_VAR:VELZ_VAR) = Ux(VELX_VAR:VELZ_VAR)-du*as*cs*sbn*Bt/d
    Ux(VELN_VAR)     = Ux(VELN_VAR)+du*af*cf/d
    Ux(PRES_VAR)     = Ux(PRES_VAR)+du*af*c2
    Ux(MAGX_VAR:MAGZ_VAR) = Ux(MAGX_VAR:MAGZ_VAR)+du*as*c*Bt/sqd

    ! Compute derivative corrections and correct states
    counter = 10 ! Maximum number of gradient reductions
    do
       Ut(DENS_VAR)     = U(VELN_VAR,i)*Ux(DENS_VAR)+U(DENS_VAR,i)*Ux(VELN_VAR)

       Ut(VELX_VAR:VELZ_VAR) = U(VELN_VAR,i)*Ux(VELX_VAR:VELZ_VAR)-U(MAGN_VAR,i)*Ux(MAGX_VAR:MAGZ_VAR)/U(DENS_VAR,i)

       Ut(VELN_VAR)     = Ut(VELN_VAR)-grav(i)*dx(i)+ &
            (Ux(PRES_VAR)+dot_product(U(MAGX_VAR:MAGZ_VAR,i),Ux(MAGX_VAR:MAGZ_VAR)))/U(DENS_VAR,i)

       Ut(PRES_VAR)     = (U(GAMC_VAR,i)*U(PRES_VAR,i))*Ux(VELN_VAR)+U(VELN_VAR,i)*Ux(PRES_VAR)

       Ut(MAGX_VAR:MAGZ_VAR) = U(MAGX_VAR:MAGZ_VAR,i)*Ux(VELN_VAR)- &
            U(MAGN_VAR,i)*Ux(VELX_VAR:VELZ_VAR)+U(VELN_VAR,i)*Ux(MAGX_VAR:MAGZ_VAR)

       Uc(:,i) = U(:,i)-0.5*dt*Ut/dx(i)
       Um(:,i) = Uc(:,i)-0.5*Ux
       Up(:,i) = Uc(:,i)+0.5*Ux


       ! We compute time derivatives of species abundances only
       ! to supply them to the equation of state in this loop,
       ! so that the internal energy is computed using consistent
       ! physical states. The actual advancement of species is
       ! done by a completely different Lagrangian algorithm.

       do ins = SPECIES_BEGIN,SPECIES_END
          Ut(ins) = 0.5*dt*U(VELN_VAR,i)*Ux(ins)/dx(i)
        
          specIndex = ins - NPROP_VARS
          Xm(specIndex) = Um(ins,i)-Ut(ins)
          Xp(specIndex) = Up(ins,i)-Ut(ins)
       end do

       ! Prevent unphysical states
       d = hy_eswitch*min(U(DENS_VAR,i-1),U(DENS_VAR,i),U(DENS_VAR,i+1))
       s = hy_eswitch*min(U(PRES_VAR,i-1),U(PRES_VAR,i),U(PRES_VAR,i+1))
       if( Um(DENS_VAR,i) > d .AND. Up(DENS_VAR,i) > d .AND. &
           Um(PRES_VAR,i) > s .AND. Up(PRES_VAR,i) > s ) then
          exit
       else
          if( counter < 0 ) then
             Ux = 0.
             Um(:,i) = Uc(:,i)
             Up(:,i) = Uc(:,i)
!!$             if (hy_meshMe.EQ.MASTER_PE) &
!!$                  print *, "[hy_8wv_interpolate] Use First order interpolation to prevent unphysical states"
             exit
          else
             Ux = 0.5*Ux
             counter = counter-1
          end if
       end if
    end do

    ! Fix abundances if needed
    xmsum = 0.
    xpsum = 0.

    do ins = SPECIES_BEGIN,SPECIES_END
       specIndex = ins - NPROP_VARS
       Xm(specIndex) = max(0.,min(1.,Xm(specIndex)))
       Xp(specIndex) = max(0.,min(1.,Xp(specIndex)))
       xmsum = xmsum+Xm(specIndex)
       xpsum = xpsum+Xp(specIndex)
    end do

    do ins = SPECIES_BEGIN,SPECIES_END
       specIndex = ins - NPROP_VARS
       Xm(specIndex) = Xm(specIndex)/xmsum
       Xp(specIndex) = Xp(specIndex)/xpsum
    end do



    ! Eos call takes inputs of density and pressure,
    ! and outputs internal energy and gamc
    vecLen=1
    eosMf =1.

    eosData(EOS_DENS)=Um(DENS_VAR,i)
    eosData(EOS_PRES)=Um(PRES_VAR,i)
    eosData(EOS_TEMP)=Uc(TEMP_VAR,i)
    call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)

    Um(EINT_VAR,i)=eosData(EOS_EINT)
    Um(TEMP_VAR,i)=eosData(EOS_TEMP)
    Um(GAMC_VAR,i)=eosData(EOS_GAMC)


    eosData(EOS_DENS)=Up(DENS_VAR,i)
    eosData(EOS_PRES)=Up(PRES_VAR,i)
    eosData(EOS_TEMP)=Uc(TEMP_VAR,i)

    call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
    Up(EINT_VAR,i)=eosData(EOS_EINT)
    Up(TEMP_VAR,i)=eosData(EOS_TEMP)
    Up(GAMC_VAR,i)=eosData(EOS_GAMC)

  end do

end subroutine hy_8wv_interpolate
