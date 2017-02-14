!!****if* source/physics/Eos/EosMain/multiTemp/Multigamma/eos_byTemp
!!
!! NAME
!!
!!  eos_byTemp
!!
!!
!! SYNOPSIS
!!
!!  call eos_byTemp ( integer(IN) :: eos_jlo,
!!                    integer(IN) :: eos_jhi,
!!          optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM),
!!    optional,target,integer(IN) :: componentMask(N_EOS_TEMP),
!!          optional,    real(IN) :: eleFrac,
!!          optional,    real(IN) :: ggProdEle )
!!
!!
!! DESCRIPTION
!!
!!  Multitemp Gamma Eos routine with extended functionality.
!!
!!  Eos computations are done for up to N_EOS_TEMP possible fluid components,
!!  which each may have a different temperature (and pressure, internal energy,
!!  etc.).
!!
!!  Given a temperature temp [K], density den [g/cm**3], and a composition 
!!  characterized by abar and zbar, this routine returns some of the other 
!!  thermodynamic quantities. Of prime interest is the pressure [erg/cm**3], 
!!  specific thermal energy [erg/gr], and their derivatives with respect
!!  to temperature and density.
!!
!!  
!!  For opacity purposes, the number density of electrons,
!!  the chemical potential of the electrons, and the pressure due
!!  to electrons is also returned. DEV: really?
!!  
!!  This routine uses planckian photons, an ideal gas of ions,
!!  and an ideal gas of electrons.
!!  
!!
!! ARGUMENTS
!!
!!  eos_jlo  -- low index of calculation range
!!  eos_jhi  -- high index of calculation range
!!  mask     --  Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of basic variables.
!!  eleFrac   - unused
!!  ggProdEle - should be 1.0/(eos_gammaEle-1.0) * eos_gasConstant
!!
!!  componentMask : selects fluid ccomponents to include in the computation.
!!                  Each mask entry should normally be either 0 (do not include)
!!                  or 1 (do include this fluid component).
!!                  Initially intended use:
!!                  componentMask(0)  determines whether to include ions
!!                  componentMask(1)  determines whether to include electrons
!!                  componentMask(2)  determines whether to include radiation
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  equation of state communication through eos_vecData
!!
!!  btemp    = temperature
!!  den      = density
!!  abar     = average number of nucleons per nuclei
!!  zbar     = average number of protons per nuclei
!!  z2bar    = square of zbar
!!  ytot1    = total number of moles per gram
!!  ye       = electron mole number
!!
!!  Since this subroutine has optional arguments, an explicit interface is needed.
!!  Calling program units should therefore have a line like
!!    use eos_mtInterface, ONLY : eos_byTemp
!!  no matter whether they call eos_byTemp with or without the optional mask argument.
!!***

subroutine eos_byTemp(eos_jlo,eos_jhi,mask,componentMask,eleFrac,ggProdEle)

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use eos_vecData, ONLY: tempRow , denRow, abarRow, zbarRow, &
       ptotRow, etotRow, stotRow, &
       dpdRow, dptRow, dstRow, dedRow, detRow, dsdRow, &
       deaRow, dezRow, & !Calhoun
       pelRow, neRow, etaRow, gamcRow, cvRow, cpRow, &
       tempRadRow, tempIonRow, tempEleRow, &
       eCompRow, pCompRow

  !! physical constants to high precision
  use eos_helmConstData, ONLY: eos_kerg, eos_kergavo, eos_ao3, eos_avo, eos_avoInv, eos_sioncon
  use Eos_data, ONLY: eos_gasConstant, eos_smallT

  implicit none


#include "constants.h" 
#include "Flash.h"  
#include "Eos.h"
#include "Eos_components.h"

!! Arguments
  integer, intent(IN) :: eos_jlo, eos_jhi
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
  integer,optional,target, dimension(N_EOS_TEMP),INTENT(in)::componentMask
  real,optional,INTENT(in) :: eleFrac,ggProdEle

!! Local variables

  character*90 ::  internalFile
  integer      ::  i,j

  real         ::  btemp,den,abar,zbar,z2bar,ytot1,ye

  real         ::  x1,x2,x3,x4,x5,x6,x7, &
                   y0,y1,y2,y3,y4, &  ! don't reuse variables -- it's confusing!
       deni,tempi,kt, & 
       prad,dpraddd,dpraddt,erad,deraddd,deraddt, &
       srad, dsraddd, dsraddt, & 
       xni,pion,dpiondd,dpiondt,eion,deiondd, deiondt,& 
       sion,dsiondd, dsiondt, & 
       pele,dpepdd,dpepdt,eele,deepdd,deepdt, &
       sele, dsepdd,dsepdt, & 
       pres,dpresdd,dpresdt,ener,denerdt, &
       entr, & 
       presi,chit,chid,gamc,kavoy

 real :: btempRad,tempiRad
 real :: btempIon,tempiIon, ktIon
 real :: btempEle,tempiEle
 real :: cv, cp, etaele, xnefer,          &
         denerdd,dentrdd,dentrdt


!!  for the interpolations
  real             xt,xd,mxt,mxd, & 
       si0t,si1t,si2t,si0mt,si1mt,si2mt, & 
       si0d,si1d,si2d,si0md,si1md,si2md, & 
       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, & 
       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, & 
       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, & 
       zFunc,z0,z1,z2,z3,z4,z5,z6, & ! Split up confusing calculations
       psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, & 
       dpsi2,ddpsi2,din,h5, & 
       xpsi0,xpsi1,h3dpd, & 
       w0t,w1t,w2t,w0mt,w1mt,w2mt, & 
       w0d,w1d,w2d,w0md,w1md,w2md
  real :: h3e, h3x
  real :: fi(36)


!!  Want NO coulomb corrections
  real ktinv,dxnidd,dsdd,lami,inv_lami,lamidd, & 
       s0,s1,s2,s3,s4, &   ! all temporary variables, reuse is confusing
       plasg,plasg_inv,plasgdd,plasgdt,a1,b1,c1,d1cc,e1cc,a2,b2,c2, & 
       ecoul,decouldd,decouldt,pcoul,dpcouldd,dpcouldt, & 
       scoul,dscouldd,dscouldt

!Added by Calhoun for calculations for the Aprox13t network
   real    :: deradda,dxnida,dpionda,deionda,dsepda,deepda,&
     &     dsda,dsdda,lamida,plasgda,denerda,deraddz,deiondz,deepdz, &
     &     dsepdz,plasgdz,denerdz
   logical :: bAprox13t  ! becomes true if variables for Aprox13t network are set


   integer,POINTER :: cMask(:)
   integer,save,target :: defaultComponentMask(N_EOS_TEMP)
   data defaultComponentMask / N_EOS_TEMP * 1 /

  real,   parameter :: third = 1.0e0/3.0e0, & 
       forth = 4.0e0/3.0e0, & 
       qe    = 4.8032068e-10,   & 
       esqu  = qe * qe


!!  for the uniform background coulomb correction
!!$  data a1,b1,c1,d1cc,e1cc,a2,b2,c2 & 
!!$       /-0.898004e0, 0.96786e0, 0.220703e0, -0.86097e0, & 
!!$       2.5269e0  , 0.29561e0, 1.9885e0,    0.288675e0/


! --------------------------------------------------------------------------

  !! ***********NO statement function declarations **********

! -----------------------------------------------------------------------------

!!  popular format statements
03 format(1x,4(a,1pe11.3))
04 format(1x,4(a,i4))

! ------------------------------------------------------------------------------

! Initial testing of masks
  bAprox13t = .false.  
  if (present(componentMask)) then
     cMask => componentMask
  else
     cMask => defaultComponentMask
  end if

!  Note that there should be things added here for entropy eventually
  bAprox13t = .false.  
  if (present(mask)) then
     if (mask(EOS_DEA).or.mask(EOS_DEZ)) &
          bAprox13t = .true.    ! we will do calculations for Aprox13t network
  end if


  !!  normal execution starts here, start of pipeline loop, no input checking
!  call Timers_start("eos_byTemp")

  do j=eos_jlo,eos_jhi

     btemp  = tempRow(j)
     den    = denRow(j)
     abar   = abarRow(j)
     zbar   = zbarRow(j)
     ytot1  = 1.0e0/abar
     ye     = ytot1 * zbar


     !!  frequent combinations
     deni    = 1.0e0/den
     tempi   = 1.0e0/btemp 
     kt      = eos_kerg * btemp
     ktinv   = 1.0e0/kt
     kavoy   = eos_kergavo * ytot1


     !!  radiation section:
     btempRad = tempRadRow(j)
     tempiRad   = 1.0e0/max(btempRad,eos_smallT)
     prad    = eos_ao3 * btempRad * btempRad * btempRad * btempRad
     dpraddt = 4.0e0 * prad * tempiRad
     dpraddd = 0.0e0

     x1      = prad * deni 
     erad    = 3.0e0 * x1
     deraddd = -erad*deni
     deraddt = 4.0e0 * erad * tempiRad
     ! Calhoun next two lines
     deradda = 0.0e0
     deraddz = 0.0e0

     srad    = (x1 + erad)*tempiRad
     dsraddd = (dpraddd*deni - x1*deni + deraddd)*tempiRad
     dsraddt = (dpraddt*deni + deraddt - srad)*tempiRad


     !!  ion section:
     btempIon = tempIonRow(j)
     tempiIon   = 1.0e0/btempIon 
     ktIon      = eos_kerg * btempIon
     dxnidd  = eos_avo * ytot1
     xni     = dxnidd * den

     pion    = xni * ktIon
     dpiondd = eos_avo * ytot1 * ktIon
     dpiondt = xni * eos_kerg
     
   

     eion    = 1.5e0 * pion * deni
     deiondd = (1.5e0 * dpiondd - eion)*deni
     deiondt = 1.5e0 * xni * eos_kerg *deni

     if (bAprox13t) then
        dxnida  = -xni*ytot1 !Calhoun
        dpionda = dxnida * ktIon !Calhoun
        deionda = 1.5e0 * dpionda*deni  !Calhoun
        deiondz = 0.0e0 !Calhoun
     end if

     !!  sackur-tetrode equation for the ion entropy of 
     !!  a single ideal gas characterized by abar
     x2      = abar*abar*sqrt(abar) * deni*eos_avoinv
     y0      = eos_sioncon * btempIon
     z0      = x2 * y0 * sqrt(y0)
     sion    = (pion*deni + eion)*tempi + kavoy*log(z0)
     dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi    &
          - kavoy * deni
     dsiondt = (dpiondt*deni + deiondt)*tempi   &
          - (pion*deni + eion) * tempi*tempi  &
          + 1.5e0 * kavoy * tempi


     !!  electron-positron section:
     btempEle = tempEleRow(j)
     tempiEle = 1.0e0/btempEle
     !!  number of electrons per volume unit is ye*den * eos_avo
     din = ye*den

     xnefer  = din * eos_avo        !DEV: Is this valid in MKS units, too? - KW

     !!  the desired electron thermodynamic quantities
     pele    = eos_gasConstant * den * btempEle * ye
     dpepdt  = eos_gasConstant * den            * ye
     dpepdd  = eos_gasConstant * bTempEle       * ye

     eele    = ggprodEle * btempEle * ye
     deepdt  = ggprodEle            * ye
     deepdd  = 0.0

     sele    = (pele * deni + eele) * tempiEle
     dsepdt  = ( (dpepdt * deni + deepdt) - (pele * deni + eele)*tempiEle ) * tempiEle
     dsepdd  = ( (dpepdd - pele * deni  ) * deni  + deepdd ) * tempiEle

     if (bAprox13t) then
        dsepda  = 0.0                                  !Calhoun ?
        dsepdz  = 0.0                                  !Calhoun ?
        deepda  = 0.0                                                     !Calhoun ?
        deepdz  = 0.0                                                 !Calhoun ?
     end if


     !!  coulomb section: GONE




     !!  sum all the components
     pres = cMask(EOSCOMP_RAD)*prad + &
          cMask(EOSCOMP_ION)*  pion + &
          cMask(EOSCOMP_ELE)* (pele)
     ener = cMask(EOSCOMP_RAD)*erad + &
          cMask(EOSCOMP_ION)*  eion + &
          cMask(EOSCOMP_ELE)* (eele)
     entr = cMask(EOSCOMP_RAD)*srad + &
          cMask(EOSCOMP_ION)*  sion + &
          cMask(EOSCOMP_ELE)* (sele)

     dpresdd = cMask(EOSCOMP_RAD)*dpraddd + cMask(EOSCOMP_ION)*dpiondd + cMask(EOSCOMP_ELE)*(dpepdd)
     dpresdt = cMask(EOSCOMP_RAD)*dpraddt + cMask(EOSCOMP_ION)*dpiondt + cMask(EOSCOMP_ELE)*(dpepdt)

     denerdd = cMask(EOSCOMP_RAD)*deraddd + cMask(EOSCOMP_ION)*deiondd + cMask(EOSCOMP_ELE)*(deepdd)
     denerdt = cMask(EOSCOMP_RAD)*deraddt + cMask(EOSCOMP_ION)*deiondt + cMask(EOSCOMP_ELE)*(deepdt)
#ifdef DEBUG_DEC2010
     print*,'denerdt,cMask(EOSCOMP_RAD)*deraddt,cMask(EOSCOMP_ION)*deiondt,cMask(EOSCOMP_ELE)*deepdt:',&
          denerdt,cMask(EOSCOMP_RAD)*deraddt,cMask(EOSCOMP_ION)*deiondt,cMask(EOSCOMP_ELE)*(deepdt)
#endif

     dentrdd = cMask(EOSCOMP_RAD)*dsraddd + cMask(EOSCOMP_ION)*dsiondd + cMask(EOSCOMP_ELE)*(dsepdd)
     dentrdt = cMask(EOSCOMP_RAD)*dsraddt + cMask(EOSCOMP_ION)*dsiondt + cMask(EOSCOMP_ELE)*(dsepdt)

     !!  form gamma_1
     presi = 1.0e0/pres
     chit  = btemp*presi * dpresdt
     chid  = dpresdd * den*presi
     x7     = pres * deni * chit/(btemp * denerdt)
     gamc  = chit*x7 + chid
     cv    = denerdt
     cp    = cv*gamc/chid



     !!  store the output -- note that many of these are not used by the calling program!
     ptotRow(j)   = pres   !used by Eos as EOS_PRES = PRES_VAR
     etotRow(j)   = ener   !used by Eos as EOS_EINT = EINT_VAR
     stotRow(j)   = entr   !this is entropy, used by Eos as EOS_ENTR

     pCompRow(EOSCOMP_ION,j) = pion
     pCompRow(EOSCOMP_ELE,j) = pele
     pCompRow(EOSCOMP_RAD,j) = prad
     eCompRow(EOSCOMP_ION,j) = eion
     eCompRow(EOSCOMP_ELE,j) = eele
     eCompRow(EOSCOMP_RAD,j) = erad

     dpdRow(j)    = dpresdd  ! used as EOS_DPD
     dptRow(j)    = dpresdt  ! used as EOS_DPT ALWAYS used by MODE_DENS_PRES in Eos.F90

     dedRow(j)    = denerdd  ! used as EOS_DED
     detRow(j)    = denerdt  ! used as EOS_DET  ALWAYS used by MODE_DENS_EI in Eos.F90

     if (bAprox13t) then
        denerda = cMask(EOSCOMP_RAD)*deradda + cMask(EOSCOMP_ION)*deionda + cMask(EOSCOMP_ELE)*(deepda)  !Cal
        denerdz = cMask(EOSCOMP_RAD)*deraddz + cMask(EOSCOMP_ION)*deiondz + cMask(EOSCOMP_ELE)*(deepdz)  !Cal
        deaRow(j)    = denerda  !Calhoun EOS_DEA
        dezRow(j)    = denerdz  !Calhoun EOS_DEZ
     end if


     dsdRow(j)    = dentrdd  ! used as EOS_DSD  
     dstRow(j)    = dentrdt  ! used as EOS_DST  

     !UNUSED     pradRow(j)   = prad
     !UNUSED     eradRow(j)   = erad
     !UNUSED     sradRow(j)   = srad

     !UNUSED     pionRow(j)   = pion
     !UNUSED     eionRow(j)   = eion
     !UNUSED     sionRow(j)   = sion
 
     pelRow(j)   = pele     ! used as EOS_PEL
     !UNUSED     eeleRow(j)   = eele
     !UNUSED     seleRow(j)   = sele

     neRow(j)    = xnefer   ! used as EOS_NE  
     !UNUSED HERE etaRow(j) = etaele     ! used as EOS_ETA 

     gamcRow(j)   = gamc    !used as EOS_GAMC = GAMC_VAR

     cvRow(j)     = cv      ! EOS_CV
     cpRow(j)     = cp      ! EOS_CP

     !!  end of vectorization loop
  enddo

!   call Timers_stop("eos_byTemp")
#ifdef DEBUG_DEC2010
  print*,'eos_byTemp ret:detRow:',detRow(eos_jlo:eos_jhi)
#endif

  return

end subroutine eos_byTemp
