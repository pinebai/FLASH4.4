!!****if* source/Grid/GridBoundaryConditions/Flash2HSE/gr_applyFlash2HSEBC
!!
!! NAME
!!
!!  gr_applyFlash2HSEBC
!!
!! SYNOPSIS
!!
!!  call gr_applyFlash2HSEBC(integer(IN) :: bcType,
!!                           integer(IN) :: bcDir,
!!                           integer(IN) :: guard,
!!                           real(INOUT) :: dataRow(2*guard,NUNK_VARS),
!!                           integer(IN) :: face,
!!                           real(IN)    :: cellCenterSweepCoord(*),
!!                           real(IN)    :: secondCoord,
!!                           real(IN)    :: thirdCoord)
!!
!! DESCRIPTION
!!
!!   Implementation of hydrostatic boundary conditions as they were found in
!!   the default FLASH2 implementation.
!!
!! ARGUMENTS
!!
!!   bcType : the type of boundary condition being applied.
!!            Should be a negative number selecting one of the variants
!!            of HSE BCs.
!!
!!   bcDir :  the dimension along which to apply boundary conditions,
!!            can take values of IAXIS, JAXIS, and KAXIS.
!!
!!   guard :  number of guardcells
!!
!!   dataRow: storage for the data being operated upon.
!!            one half is input, the other half is output.
!!            Which is which depends on face.
!!
!!   face :   can take values LOW and HIGH, defined in constants.h,
!!            to indicate whether to apply boundary on lowerface or 
!!            upperface
!!
!!  cellCenterSweepCoord : vector of (at least) 2*guard cell center coordinate
!!                         values in the bcDir direction. The elements of this
!!                         array give the locations of corresponding data in
!!                         dataRow.
!!
!!  secondCoord,thirdCoord: scalar cell center coordinate values in the coordinate
!!                         directions perpendicular to the direction given by bcDir.
!!                         This is used in this implementation only for passing
!!                         the information along to Gravity_accelAtCoords, which
!!                         is the Gravity interface that is invoked to get the
!!                         strength of the gravitational field.
!!                         The meaning depends on the sweep direction bcDir as
!!                         follows
!!                          bcDir   |    secondCoord       thirdCoord
!!                          ------------------------------------------
!!                          IAXIS   |    Y(j) *            Z(k) **
!!                          JAXIS   |    X(i)              Z(k) **
!!                          KAXIS   |    X(i)              Y(j)
!!                         *)  if NDIM > 1
!!                         **) if NDIM > 2
!!                         These dummy arguments are ignored (and an implementation
!!                         of this interface should not attempt to access them) if
!!                         they do not make sense based on the dimensionality of
!!                         the problem.
!!
!! NOTES
!!
!!
!!  This implementation may or may not be appropriate for a given application.
!!  Most setups in FLASH2 that used hydrostatic BCs actually implemented their
!!  own version instead of using the default implementation from which this
!!  file is derived.
!!
!!***

subroutine gr_applyFlash2HSEBC(bcType,bcDir,guard,dataRow,face,&
     cellCenterSweepCoord, secondCoord,thirdCoord)

  use Driver_interface, ONLY : Driver_abortFlash
  use Gravity_interface, ONLY : Gravity_accelAtCoords
  use Eos_interface, ONLY : Eos
  use ut_interpolationInterface, ONLY : ut_quadraticInterpol,&
       ut_quadraticCellAverageInterpol
  use Grid_data, ONLY : gr_smalle, gr_iguard, gr_jguard, gr_kguard
  use gr_bcData, ONLY : gr_bcEintSwitch,gr_bcUseGravity

  implicit none

#include "Flash.h"
#include "constants.h"      
#include "Eos.h"

  integer,intent(IN):: bcType,bcDir,guard,face
  real,dimension(2*guard,NUNK_VARS),intent(INOUT)::dataRow
  real,intent(IN):: cellCenterSweepCoord(*), secondCoord,thirdCoord

  integer :: ic, ib, ie

  integer :: q

  real, dimension(max(2*gr_iguard, 2*gr_jguard*K2D, 2*gr_kguard*K3D))    :: r, p, xzn, xznr, xznl, tg, t, ei, rho
  real :: pc, pl, pr, m
  real :: dxzn, dxtoi, dxip1, xcl, xcr
  real :: rint, pint, gint, xint, rgc, rgcl, rgcr, rgl, rgr, pcl, pcr
  real :: go, ro
!!unused?!  real :: abar,zbar,dpt,dpd,det,ded,kgc,pel,ne,eta
  real :: otherCoord1(1), otherCoord2(1)

  real, dimension(guard):: energyKinetic,energyInternal
  real, dimension(NSPECIES*guard) :: massFraction
  real, dimension(EOS_NUM*guard) :: eosData

  integer :: i, eosVecLen,pres,dens,gamc,temp,abar,zbar,eint
  integer :: normalVelocityVar

!    controls the interpolation for the HSE BCs
!
  integer, parameter :: CONST = 0, LINEAR = 1, QUADRATIC = 2, INVSQ = -2
  integer, parameter :: interp = QUADRATIC
!
!    controls how bnd_hydrostat handles velocities
!
  integer, parameter :: VEL_REFLECT = 0, VEL_OUTFLOW = 1, VEL_DIODE = 2
  integer, save :: hsevel = VEL_REFLECT
  

#ifndef FLASH_EOS
  call Driver_abortFlash('Cannot execute gr_applyFlash2HSEBC without the Eos unit!')
#endif

  q = max(2*gr_iguard, 2*gr_jguard*K2D, 2*gr_kguard*K3D)

  if (bcType.eq.HYDROSTATIC_F2_NVDIODE .OR. bcType.eq.HYDROSTATIC_NVDIODE) then
     hsevel = VEL_DIODE
  else if (bcType.eq.HYDROSTATIC_F2_NVOUT .OR. bcType.eq.HYDROSTATIC_NVOUT) then
     hsevel = VEL_OUTFLOW
  else 
     hsevel = VEL_REFLECT
  endif

  if (bcDir==IAXIS) then
     normalVelocityVar = VELX_VAR
  else if (bcDir==JAXIS) then
     normalVelocityVar = VELY_VAR
  else
     normalVelocityVar = VELZ_VAR
  end if

  if (gr_bcUseGravity) then
     otherCoord1(1) = secondCoord
     otherCoord2(1) = thirdCoord
  end if

!               -{X,Y,Z} hydrostatic boundary

  if (face==LOW) then

     ic = guard+1

     xzn(1:2*guard) = cellCenterSweepCoord(1:2*guard)
     xznl(1) = xzn(1) - 0.5*(xzn(2) - xzn(1))
     xznl(2:2*guard) = 0.5 * (xzn(1:2*guard-1) + xzn(2:2*guard))
     xznr(1:2*guard-1) = xznl(2:2*guard)
     xznr(2*guard) = xzn(2*guard) + 0.5*(xzn(2*guard) - xzn(2*guard - 1))
     if (gr_bcUseGravity) then
        if (bcDir==IAXIS) then
           call Gravity_accelAtCoords(2*guard, xzn, otherCoord1(1:1), otherCoord2(1:1), bcDir, tg)
        else if (bcDir==JAXIS) then
           call Gravity_accelAtCoords(2*guard, otherCoord1(1:1), xzn, otherCoord2(1:1), bcDir, tg)
        else
           call Gravity_accelAtCoords(2*guard, otherCoord1(1:1), otherCoord2(1:1), xzn, bcDir, tg)
        end if
     else
        tg = 0.0
     endif

     p    = dataRow(1:2*guard,PRES_VAR)
     r    = dataRow(1:2*guard,DENS_VAR)
     t    = dataRow(1:2*guard,TEMP_VAR)
!!!     ei   = dataRow(1:2*guard,EINT_VAR)

!
! find p, rho, g at the interface to the first boundary cell
!
     xint = xznl(ic)

     call lagrange_cc_interp(pint,xint, p(ic),xzn(ic),xznl(ic))
     call lagrange_cc_interp(rint,xint, r(ic),xzn(ic),xznl(ic))
     call lagrange_interp   (gint,xint,tg(ic),xzn(ic))

     if (interp .eq. INVSQ) then 
        if (tg(ic) .eq. tg(ic+1)) then 
           print *,'error -- using const grav and INVSQ in gr_applyFlash2HSEBC'
           call Driver_abortFlash('should not be using INVSQ in gr_applyFlash2HSEBC')
        else
           go = (xzn(ic)-xzn(ic+1))/(1./sqrt(-tg(ic)) - 1./sqrt(-tg(ic+1)))
           go = -go*go
           ro = xzn(ic) - sqrt(go/tg(ic))
           gint = go/((xint-ro)*(xint-ro))
        endif
     endif


     do i = guard, 1, -1

        dataRow(i,:)     = dataRow(ic,:)

        select case (hsevel) 
        case (VEL_REFLECT)
           dataRow(i,normalVelocityVar) = -dataRow(guard+ic-i,normalVelocityVar)
        case (VEL_OUTFLOW)
           dataRow(i,normalVelocityVar) =  dataRow(ic,normalVelocityVar)
        case (VEL_DIODE)
           dataRow(i,normalVelocityVar) =  min(dataRow(ic,normalVelocityVar),0.)
        end select

        r(i)    = rint
        t(i)    = dataRow(ic,TEMP_VAR)
!!        xn(i,:) = dataRow(i,SPECIES_BEGIN:SPECIES_END)

        dxzn = xznr(i) - xznl(i)

        if (interp .eq. CONST) then
           pr = pint
           pc = pint - rint*gint * dxzn*0.5
           pl = pint - rint*gint * dxzn

           p(i) = pc
        else if (interp .eq. LINEAR) then
           dxip1 = xzn (i+2) - xzn (i+1)

           m = (tg(i+2) - tg(i+1))/dxip1

           pr = pint
           pc = pint-rint*(gint*dxzn*0.5 - 0.125*m*dxzn*dxzn)
           pl = pint-rint*(gint*dxzn     - 0.5  *m*dxzn*dxzn)

           p(i) = (pr + 4.*pc + pl)/6.
        else if (interp .eq. QUADRATIC) then
           xcl = xint -    dxzn/3.
           xcr = xint - 2.*dxzn/3.

           rgr = 0.
           ! It is not clear why the following extrapolation should be used
           ! as it is, why not use interpolation where possible?
           ! e.g., could use {tg,xzn}(max(i-1, 2))     - KW
           call lagrange_interp_integrate(rgcr,xcr,xint,tg(i+1),xzn(i+1))
           call lagrange_interp_integrate(rgcl,xcl,xint,tg(i+1),xzn(i+1))
           call lagrange_interp_integrate(rgl,xznl(I),xint,tg(i+1),xzn(i+1))

           pr  = pint - rint*rgr
           pcl = pint - rint*rgcl
           pcr = pint - rint*rgcr
           pl  = pint - rint*rgl

           p(i) = (pr + 3.*pcl + 3.*pcr + pl)/8.
        else if (interp .eq. INVSQ) then
           p(i) = pint + rint*go/(xint-ro)    + rint/dxzn*go*log((xznl(i)-ro)/(xint-ro))
           pl   = pint - rint*go/(xznl(i)-ro) + rint*go/(xint-ro)
        endif

        pint = pl
        rint = rint
        xint = xznl(i)

        if (interp .eq. CONST) then
           gint = tg(max(i-1,1))
        else if (interp .eq. LINEAR) then
           m = (tg(i+1) - tg(i))/(xzn(i+1)-xzn(i))
           gint  = tg(i) + m*(xint - xzn(i))
        else if (interp .eq. QUADRATIC) then
           call lagrange_interp(gint,xint,tg(i),xzn(i))
        else if (interp .eq. INVSQ) then
           gint = go/((xint-ro)*(xint-ro))
        endif

     enddo

     ib = 1
     ie = guard
     eosVecLen = guard
     ! These integers are indices into the location in eosData just before the storage area for the appropriate variable.
     pres = (EOS_PRES-1)*eosVecLen
     dens = (EOS_DENS-1)*eosVecLen
     temp = (EOS_TEMP-1)*eosVecLen
     gamc = (EOS_GAMC-1)*eosVecLen
     eint = (EOS_EINT-1)*eosVecLen
     abar = (EOS_ABAR-1)*eosVecLen
     zbar = (EOS_ZBAR-1)*eosVecLen

     !! Fill up two scratch arrays. 
     !! energyKinetic holds velocity vector information -- 1/2 * Vmag**2
     !! energyInternal holds eint (directly)  or energyTotal - ekinetic (calculated),
     !!          depending upon eintSwitch

     energyKinetic(1:eosVecLen) = 0.5*&
          (dataRow(ib:ie,VELX_VAR)**2 +  & 
          dataRow(ib:ie,VELY_VAR)**2 +  & 
          dataRow(ib:ie,VELZ_VAR)**2)
     energyInternal(1:eosVecLen) = dataRow(ib:ie,EINT_VAR)

     do i = 1,eosVecLen
        if (dataRow(ib+i-1,ENER_VAR) > &
             (1.+ gr_bcEintSwitch)*energyKinetic(i)) then
           energyInternal(i) = dataRow(ib+i-1,ENER_VAR) - energyKinetic(i)
        end if
        energyInternal(i) = max(energyInternal(i), gr_smalle)
        massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
             dataRow(ib+i-1,SPECIES_BEGIN:SPECIES_END)
     end do

     eosData(pres+1:pres+eosVecLen) =       p(ib:ie)
     eosData(dens+1:dens+eosVecLen) =       r(ib:ie)
     eosData(temp+1:temp+eosVecLen) =       t(ib:ie)
     eosData(gamc+1:gamc+eosVecLen) = dataRow(ib:ie,GAMC_VAR)
     eosData(eint+1:eint+eosVecLen) = energyInternal(1:eosVecLen)

     !!  These values may be overwritten within the internal Eos routines.
     !!  They are added here for cases where the mass fractions will not be used
     !!  to calculate abar and zbar.
#ifdef USE_EOS_YE

#ifdef YE_MSCALAR
     !! cal says abar=1/sumy
     !! cal says zbar=ye / sumy and he claims sumy are never zero
     eosData(abar+1:abar+eosVecLen) =  1.0 / dataRow(ib:ie,SUMY_MSCALAR)
     eosData(zbar+1:zbar+eosVecLen) = dataRow(ib:ie,YE_MSCALAR) / dataRow(ib:ie,SUMY_MSCALAR)
#else
     call Driver_abortFlash("gr_applyFlash2HSEBC called in USE_EOS_YE mode, but no YE_MSCALAR is defined")
#endif

#endif


     call Eos(MODE_DENS_PRES, guard, eosData, massFraction)
     do i=1,guard
        dataRow(i,DENS_VAR) = eosData(dens+i) ! r(i)
        dataRow(i,PRES_VAR) = eosData(pres+i) ! p(i)
        dataRow(i,TEMP_VAR) = eosData(temp+i) ! t(i) 
        dataRow(i,EINT_VAR) = eosData(eint+i) ! ei(i)
        dataRow(i,ENER_VAR) =dataRow(i,EINT_VAR) + 0.5 * (               &
             &                dataRow(i,VELX_VAR)*dataRow(i,VELX_VAR)+      &
             &                dataRow(i,VELY_VAR)*dataRow(i,VELY_VAR)+      &
             &                dataRow(i,VELZ_VAR)*dataRow(i,VELZ_VAR))
     enddo

  else if (face==HIGH) then

!               +{X,Y,Z} hydrostatic boundary
     
     ic = guard

     xzn(1:2*guard) = cellCenterSweepCoord(1:2*guard)
     xznl(1) = xzn(1) - 0.5*(xzn(2) - xzn(1))
     xznl(2:2*guard) = 0.5 * (xzn(1:2*guard-1) + xzn(2:2*guard))
     xznr(1:2*guard-1) = xznl(2:2*guard)
     xznr(2*guard) = xzn(2*guard) + 0.5*(xzn(2*guard) - xzn(2*guard - 1))
     if (gr_bcUseGravity) then
        if (bcDir==IAXIS) then
           call Gravity_accelAtCoords(2*guard, xzn, otherCoord1(1:1), otherCoord2(1:1), bcDir, tg)
        else if (bcDir==JAXIS) then
           call Gravity_accelAtCoords(2*guard, otherCoord1(1:1), xzn, otherCoord2(1:1), bcDir, tg)
        else
           call Gravity_accelAtCoords(2*guard, otherCoord1(1:1), otherCoord2(1:1), xzn, bcDir, tg)
        end if
     else
        tg = 0.0
     endif

     p    = dataRow(:,PRES_VAR)
     r    = dataRow(:,DENS_VAR)
     t    = dataRow(:,TEMP_VAR)
!!!     ei   = dataRow(:,EINT_VAR)
     !
     ! find p, rho, g at the interface to the first boundary cell
     !
     xint = xznr(ic)

     call lagrange_cc_interp(pint,xint, p(ic-2),xzn(ic-2),xznl(ic-2))
     call lagrange_cc_interp(rint,xint, r(ic-2),xzn(ic-2),xznl(ic-2))
     call lagrange_interp   (gint,xint,tg(ic-2),xzn(ic-2))

     if (interp .eq. INVSQ) then 
        if (tg(ic) .eq. tg(ic-1)) then 
           print *,'error -- using const grav and INVSQ in gr_applyFlash2HSEBC'
           call Driver_abortFlash('should not be using INVSQ in gr_applyFlash2HSEBC')
        else
           go = (xzn(ic)-xzn(ic-1))/(1./sqrt(-tg(ic)) - 1./sqrt(-tg(ic-1)))
           go = -go*go
           ro = xzn(ic) - sqrt(go/tg(ic))
           gint = go/((xint-ro)*(xint-ro))
        endif
     endif

     do i = guard+1, 2*guard

        dataRow(i,:) = dataRow(ic,:)
        select case (hsevel) 
        case (VEL_REFLECT)
           dataRow(i,normalVelocityVar) = -dataRow(2*ic+1-i,normalVelocityVar)
        case (VEL_OUTFLOW)
           dataRow(i,normalVelocityVar) =  dataRow(ic,normalVelocityVar)
        case (VEL_DIODE)
           dataRow(i,normalVelocityVar) =  max(dataRow(ic,normalVelocityVar),0.)
        end select

!!        xn(i,:) = dataRow(i,SPECIES_BEGIN:SPECIES_END)
        r(i)    = rint
        t(i)    = dataRow(ic,TEMP_VAR)

        dxzn = (xznr(i) - xznl(i))


        if (interp .eq. CONST) then
           pl = pint
           pc = pint + rint*gint * dxzn*0.5
           pr = pint + rint*gint * dxzn

           p(i) = pc
        else if (interp .eq. LINEAR) then
           dxip1 = xzn (i-2) - xzn (i-1)

           m = (tg(i-2) - tg(i-1))/dxip1

           pl = pint
           pc = pint+rint*(gint*dxzn*0.5 + 0.125*m*dxzn*dxzn)
           pr = pint+rint*(gint*dxzn     + 0.5  *m*dxzn*dxzn)

           p(i) = (pr+4.*pc+pl)/6.
        else if (interp .eq. QUADRATIC) then
           xcl = xint + 1.*dxzn/3.
           xcr = xint + 2.*dxzn/3.

           rgl = 0.
           ! Changed to (I-3) from (i-2) in the next three statements,
           ! basically for consistency with the -X direction.
           ! F2 tot_bnd handled different directions differently, probably
           ! by mistake.
           ! It is not clear why the following extrapolation should be used
           ! as it is, why not use interpolation where possible?
           ! e.g., could use {tg,xzn}(min(i-1, 2*guard-3))     - KW
           call lagrange_interp_integrate(rgcl,xcl,xint,tg(I-3),xzn(I-3))
           call lagrange_interp_integrate(rgcr,xcr,xint,tg(I-3),xzn(I-3))
           call lagrange_interp_integrate(rgr,xznr(i),xint,tg(I-3),xzn(I-3))

           pr  = pint - rint*rgr
           pcl = pint - rint*rgcl
           pcr = pint - rint*rgcr
           pl  = pint - rint*rgl

           p(i) = (pr + 3.*pcl + 3.*pcr + pl)/8.
        else if (interp .eq. INVSQ) then
           p(i) = pint + rint*go/(xint-ro) - rint/dxzn*go*log((xznr(i)-ro)/(xint-ro))
           pr   = pint + rint*go/(xint-ro) - rint*go/(xznr(i)-ro)
        endif

        pint = pr
        rint = rint
        xint = xznr(i)

        if (interp .eq. CONST) then
           gint = tg(min(i+1,guard*2))
        else if (interp .eq. LINEAR) then
           m = (tg(i-1) - tg(i))/(xzn(i-1)-xzn(i))
           gint  = tg(i) + m*(xint - xzn(i))
        else if (interp .eq. QUADRATIC) then
           call lagrange_interp(gint,xint,tg(i-2),xzn(i-2))
        else if (interp .eq. INVSQ) then
           gint = go/((xint-ro)*(xint-ro))
        endif

     enddo

     ib = guard + 1
     ie = 2*guard
     eosVecLen = guard
     ! These integers are indices into the location in eosData just before the storage area for the appropriate variable.
     pres = (EOS_PRES-1)*eosVecLen
     dens = (EOS_DENS-1)*eosVecLen
     temp = (EOS_TEMP-1)*eosVecLen
     gamc = (EOS_GAMC-1)*eosVecLen
     eint = (EOS_EINT-1)*eosVecLen
     abar = (EOS_ABAR-1)*eosVecLen
     zbar = (EOS_ZBAR-1)*eosVecLen

     !! Fill up two scratch arrays. 
     !! energyKinetic holds velocity vector information -- 1/2 * Vmag**2
     !! energyInternal holds eint (directly)  or energyTotal - ekinetic (calculated),
     !!          depending upon eintSwitch

     energyKinetic(1:eosVecLen) = 0.5*&
          (dataRow(ib:ie,VELX_VAR)**2 +  & 
          dataRow(ib:ie,VELY_VAR)**2 +  & 
          dataRow(ib:ie,VELZ_VAR)**2)
     energyInternal(1:eosVecLen) = dataRow(ib:ie,EINT_VAR)

     do i = 1,eosVecLen
        if (dataRow(ib+i-1,ENER_VAR) > &
             (1.+ gr_bcEintSwitch)*energyKinetic(i)) then
           energyInternal(i) = dataRow(ib+i-1,ENER_VAR) - energyKinetic(i)
        end if
        energyInternal(i) = max(energyInternal(i), gr_smalle)
        massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
             dataRow(ib+i-1,SPECIES_BEGIN:SPECIES_END)
     end do

     eosData(pres+1:pres+eosVecLen) =       p(ib:ie)
     eosData(dens+1:dens+eosVecLen) =       r(ib:ie)
     eosData(temp+1:temp+eosVecLen) =       t(ib:ie)
     eosData(gamc+1:gamc+eosVecLen) = dataRow(ib:ie,GAMC_VAR)
     eosData(eint+1:eint+eosVecLen) = energyInternal(1:eosVecLen)

     !!  These values may be overwritten within the internal Eos routines.
     !!  They are added here for cases where the mass fractions will not be used
     !!  to calculate abar and zbar.
#ifdef USE_EOS_YE

#ifdef YE_MSCALAR
     !! cal says abar=1/sumy
     !! cal says zbar=ye / sumy and he claims sumy are never zero
     eosData(abar+1:abar+eosVecLen) =  1.0 / dataRow(ib:ie,SUMY_MSCALAR)
     eosData(zbar+1:zbar+eosVecLen) = dataRow(ib:ie,YE_MSCALAR) / dataRow(ib:ie,SUMY_MSCALAR)
#else
     call Driver_abortFlash("gr_applyFlash2HSEBC called in USE_EOS_YE mode, but no YE_MSCALAR is defined")
#endif

#endif



     call Eos(MODE_DENS_PRES, guard, eosData, massFraction)
     do i=guard+1,2*guard
        dataRow(i,DENS_VAR) = eosData(dens+i-guard) ! r(i)
        dataRow(i,PRES_VAR) = eosData(pres+i-guard) ! p(i)
        dataRow(i,TEMP_VAR) = eosData(temp+i-guard) ! t(i)
        dataRow(i,EINT_VAR) = eosData(eint+i-guard) ! ei(i)
        dataRow(i,ENER_VAR) =dataRow(i,EINT_VAR) + 0.5 * (               &
             &                dataRow(i,VELX_VAR)*dataRow(i,VELX_VAR)+      &
             &                dataRow(i,VELY_VAR)*dataRow(i,VELY_VAR)+      &
             &                dataRow(i,VELZ_VAR)*dataRow(i,VELZ_VAR))
     enddo

  else
     call Driver_abortFlash('[gr_applyFlash2HSEBC] face is neither LOW nor HIGH!')
  endif


contains
!  
!  quadratic interpolation routines
!

  subroutine lagrange_interp(yint, xint, y, x)
    
    implicit none
    real, intent(in)    :: xint           ! where function is to be eval'd
    real, intent(in)    :: y(*), x(*)     ! function values and locations
    real, intent(out)   :: yint           ! function evaluated at xint
    
    call ut_quadraticInterpol(x, y, xint, yint)
    return
  end subroutine lagrange_interp
  
  
  subroutine lagrange_cc_interp(yint, xint, y, x, xl)
    implicit none
    real, intent(in)    :: xint           ! where function is to be eval'd
    real, intent(in)    :: y(*), x(*), xl(*) ! function values and locations
    real, intent(out)   :: yint           ! function evaluated at xint
    
    call ut_quadraticCellAverageInterpol(x, xl, y, xint, yint)
    return
  end subroutine lagrange_cc_interp
  
  
  subroutine lagrange_interp_integrate(yint, xa, xb, y, x)
    implicit none
    real, intent(in)    :: xa,xb          ! range of integration
    real, intent(in)    :: y(*), x(*)     ! function values and locations
    real, intent(out)   :: yint           ! function evaluated at xint
    
    
    real :: dx12, dx23, dx13
    real :: a, b, c
    real, parameter :: onethird = 1./3.
    
    
    dx12 = x(1)-x(2)
    dx23 = x(2)-x(3)
    dx13 = x(1)-x(3)
    
    a =     onethird*xb*xb*xb - .5*xb*xb*(x(2)+x(3)) + xb*(x(2)*x(3))
    a = a -(onethird*xa*xa*xa - .5*xa*xa*(x(2)+x(3)) + xa*(x(2)*x(3)))
    a = a / (dx12*dx13)
    
    b =     onethird*xb*xb*xb - .5*xb*xb*(x(1)+x(3)) + xb*(x(1)*x(3))
    b = b -(onethird*xa*xa*xa - .5*xa*xa*(x(1)+x(3)) + xa*(x(1)*x(3)))
    b = b / (-dx12*dx23)
    
    c =     onethird*xb*xb*xb - .5*xb*xb*(x(1)+x(2)) + xb*(x(1)*x(2))
    c = c -(onethird*xa*xa*xa - .5*xa*xa*(x(1)+x(2)) + xa*(x(1)*x(2)))
    c = c / (dx13*dx23)
    
    yint = y(1)*a + y(2)*b + y(3)*c
    return
  end   subroutine lagrange_interp_integrate
  
end subroutine gr_applyFlash2HSEBC
