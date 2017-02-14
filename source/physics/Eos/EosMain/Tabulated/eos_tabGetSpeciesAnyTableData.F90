!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabGetSpeciesAnyTableData
!!
!! NAME
!!
!!  eos_tabGetSpeciesAnyTableData
!!
!! SYNOPSIS
!!
!!  call eos_tabGetSpeciesAnyTableData(integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: maxComp,
!!                                    integer (in)  :: selTT,
!!                                    integer (in)  :: needDerivs,
!!                                    real    (out) :: resultTT(:,:))
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation an average ionization, internal energy, or heat
!!  capacity value for all components for
!!  the specified temperature, density, and species. The interpolation method is
!!  bilinear.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   maxComp :          : The highest component (ions, electrons, or matter) for which
!!                        data are to be returned
!!   selTT              : Select ZF, EN, or HC table type
!!   resultTT          : The value of the determined average ionization etc.
!!
!!***

#include "constants.h"
#include "Eos.h"

subroutine eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        wanted, &
                                        selTT,              &
                                        needDerivs,         &
                                        i,j,k,l,&
                                        tauIn, deltaIn, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        resultTT, &
                                        T1in,T2in,D1in,D2in)

  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_data,         ONLY : eos_meshMe
  use eos_tabData,      ONLY : one,ten,                    &
                               EOS_TAB_NCOMP,              &
                               EOS_TAB_NDERIVS,            &
                               EOS_TABINT_DERIV_0,         &
                               EOS_TABINT_DERIV_DT,        &
                               EOS_TABINT_DERIV_DD,        &
                               EOS_TAB_NALLTAB,            &
                               EOS_TABVT_ZF,               &
                               EOS_TABVT_EN,               &
                               EOS_TABVT_HC,               &
                               EOS_TABVT_PR,               &
                               EOS_TABVT_ENTR,             &
                                eos_useLogTables,          &
                                eosT_varTableGroupPT,      &
                                eos_allTab


  implicit none

  integer, intent (in)  :: species
  logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  integer, intent (in)  :: selTT
  integer, intent (in)  :: needDerivs
  real,    intent (out) :: resultTT(0:,:)

  type(eosT_varTableGroupPT),pointer :: thisTypeTable

  logical,intent(IN) :: lowerBoundaryDens
  logical :: upperBoundaryDens
  logical,intent(IN) :: withinBoundaryDens
  logical,intent(IN) :: lowerBoundaryTemp
  logical :: upperBoundaryTemp
  logical,intent(IN) :: withinBoundaryTemp
  logical :: found

  integer :: d,t
  integer,intent(IN) :: i,j,k,l
  integer,save :: iSave=1,kSave=1
  integer :: varType
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: g

  real :: D1,D2,T1,T2
  real,intent(IN) :: D1in,D2in,T1in,T2in
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityTT
  real :: speciesTemperatureTT
  real,intent(IN) :: tauIn,deltaIn
  real :: tau,delta

  real,parameter :: largeExponent = 77.0      ! was 0.25 * log10(HUGE(largeExponent))
  real,parameter :: largeMultiplier = 10.0**77  ! was 10.0**largeExponent
!
!
!   ...Set the handle to the ZF, EN, PR, or HC tables and associated data arrays.
!
!
  select case (selTT)
  case(EOS_TABULAR_Z)
     varType = EOS_TABVT_ZF
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_E)
     varType = EOS_TABVT_EN
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_C)
     varType = EOS_TABVT_HC
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_P)
     varType = EOS_TABVT_PR
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_S)
     varType = EOS_TABVT_ENTR
     thisTypeTable => eos_allTab(species)%tg(varType)
  case default
     varType = 0
  end select

  if (varType == 0) then
      call Driver_abortFlash ('[eos_tabGetSpeciesAnyTableData] ERROR: no handle to tables')
  end if
!
!
!   ...Get the current temperature and ion number density of the species and find:
!
!        1) the temperature/density quadrant (T1,T2,D1,D2) boundary containing
!           containing the species's temperature and density (x)
!
!        2) the table index quadrant (i,j,k,l) of the boundary, as indicated
!           in the figure below.
!
!        3) the four tabulated values (o1,o2,o3,o4)
!
!
!                    o3------------o4   T2 (j)
!                     |            |
!                     |            |
!                     |        x   |
!                     |            |
!                     |            |
!                    o1 -----------o2   T1 (i)
!
!                  D1 (k)         D2 (l)
!
!
!      In case the temperature and/or density of the species lay outside
!      the tabulated boundaries, take the corresponding boundary values
!      of the tables. The criteria as to when the species's temperature and
!      density belong within the boundary are:
!
!                         T1 =< speciesTemperature =< T2
!                         D1 =<   speciesDensity   =< D2
!
!


  T1 = T1in
  T2 = T2in
  tau = tauIn
  D1 = D1in
  D2 = D2in
  delta = deltaIn

  upperBoundaryTemp = .NOT. (lowerBoundaryTemp .OR. withinBoundaryTemp)
  upperBoundaryDens = .NOT. (lowerBoundaryDens .OR. withinBoundaryDens)

#ifdef DEBUG_EOS
  if (selTT==EOS_TABULAR_E) then
     print*,'selTT,wanted:',selTT,wanted
  end if
#endif
  do g = 1,EOS_TAB_NCOMP
     if (.NOT.wanted(g,varType)) cycle
     if (.NOT.associated(thisTypeTable%table(g)%table)) then
        if (eos_meshMe==MASTER_PE) then
           print*,'Table not associated:g,selTT',g,selTT
        end if
        cycle
     end if
     o1 = thisTypeTable%table(g)%table(i,k)
     o2 = thisTypeTable%table(g)%table(i,l)
     o3 = thisTypeTable%table(g)%table(j,k)
     o4 = thisTypeTable%table(g)%table(j,l)
     if (.NOT. withinBoundaryTemp .AND. &
          (selTT==EOS_TABULAR_E .OR. selTT==EOS_TABULAR_P) ) then
        if (eos_useLogTables) then
           speciesTemperatureTT = log10 (speciesTemperature)
        else
           speciesTemperatureTT = speciesTemperature
        end if
        if ( (lowerBoundaryTemp .AND. (eos_useLogTables .OR. T2 > 0.0)) ) then
           if (eos_useLogTables) then
              T1 = T2 - largeExponent
              tau = (speciesTemperatureTT - T1) / largeExponent
              o1 = o3 - largeExponent
              o2 = o4 - largeExponent
           else
              T1 = 0.0
              tau = speciesTemperatureTT / T2
              o1 = 0.0
              o2 = 0.0
           end if
        else if ( (upperBoundaryTemp .AND. (eos_useLogTables .OR. T2 > 0.0)) ) then
           if (eos_useLogTables) then
              T2 = T1 + largeExponent
              tau = (speciesTemperatureTT - T1) / largeExponent
              o3 = o1 + largeExponent
              o4 = o2 + largeExponent
           else
              T2 = T1 * largeMultiplier
              tau = (speciesTemperatureTT - T1) / (T2 - T1)
              o3 = o1 * largeMultiplier
              o4 = o2 * largeMultiplier
           end if
        end if
     end if
     if (.NOT. withinBoundaryDens .AND. &
          (selTT==EOS_TABULAR_P) ) then
        if (eos_useLogTables) then
           speciesDensityTT     = log10 (speciesDensity)
        else
           speciesDensityTT     = speciesDensity
        end if
        if ( (lowerBoundaryDens .AND. (eos_useLogTables .OR. D2 > 0.0)) ) then
           if (eos_useLogTables) then
              D1 = D2 - largeExponent
              delta = (speciesDensityTT - D1) / largeExponent
              o1 = o2 - largeExponent
              o3 = o4 - largeExponent
           else
              D1 = 0.0
              delta = speciesDensityTT / D2
              o1 = 0.0
              o3 = 0.0
           end if
        else if ( (upperBoundaryDens .AND. (eos_useLogTables .OR. D2 > 0.0)) ) then
           if (eos_useLogTables) then
              D2 = D1 + largeExponent
              delta = (speciesDensityTT - D1) / largeExponent
              o2 = o1 + largeExponent
              o4 = o3 + largeExponent
           else
              D2 = D1 * largeMultiplier
              delta = (speciesDensityTT - D1) / (D2 - D1)
              o2 = o1 * largeMultiplier
              o4 = o3 * largeMultiplier
           end if
        end if
     end if
     !
     !
     !   ...Do the bilinear interpolation:
     !
     !                result =   o1 * [(1-tau)*(1-delta)]
     !                         + o2 * [delta*(1-tau)]
     !                         + o3 * [tau*(1-delta)]
     !                         + o4 * [delta*tau]
     !
     !
!!$     tau   = (speciesTemperatureTT - T1) / (T2 - T1)
!!$     delta = (speciesDensityTT     - D1) / (D2 - D1)

     f1 = (one - tau) * (one - delta)
     f2 = delta * (one - tau)
     f3 = tau * (one - delta)
     f4 = delta * tau

     resultTT(EOS_TABINT_DERIV_0,g) = o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4
     !   ...Convert logarithmic form to real form (if needed).
     if (eos_useLogTables) then
        resultTT(EOS_TABINT_DERIV_0,g) = ten ** resultTT(EOS_TABINT_DERIV_0,g)
     end if

     if (needDerivs .GE. EOS_TABINT_DERIV_DT) then
        if (.NOT.withinBoundaryTemp .AND. &
            .NOT.( (selTT==EOS_TABULAR_E .OR. selTT==EOS_TABULAR_P) .AND. &
                   (  (lowerBoundaryTemp .AND. (eos_useLogTables .OR. T2 > 0.0) ) .OR. &
                      (upperBoundaryTemp .AND. (eos_useLogTables .OR. T2 > 0.0) ) ) )) then
           resultTT(EOS_TABINT_DERIV_DT,g) = 0.0
        else
           f1 =   delta - one
           f2 = - delta
           f3 =   one - delta
           f4 =   delta
           resultTT(EOS_TABINT_DERIV_DT,g) = &
                (o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4) / (T2 - T1)
           !   ...Convert logarithmic form to real form (if needed).
           if (eos_useLogTables) then
              resultTT(EOS_TABINT_DERIV_DT,g) = &
                   resultTT(EOS_TABINT_DERIV_0,g) * resultTT(EOS_TABINT_DERIV_DT,g) &
                   / speciesTemperature
           end if
        end if
     end if

     if (needDerivs .GE. EOS_TABINT_DERIV_DD) then
        if (.NOT.withinBoundaryDens .AND. &
            .NOT.( (selTT==EOS_TABULAR_P) .AND. &
                   (  (lowerBoundaryDens .AND. (eos_useLogTables .OR. D2 > 0.0) ) .OR. &
                      (upperBoundaryDens .AND. (eos_useLogTables .OR. D2 > 0.0) ) ) )) then
           resultTT(EOS_TABINT_DERIV_DD,g) = 0.0
        else
           f1 =   tau - one
           f2 =   one - tau
           f3 = - tau
           f4 =   tau
           resultTT(EOS_TABINT_DERIV_DD,g) = &
                (o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4) / (D2 - D1)
           !   ...Convert logarithmic form to real form (if needed).
           if (eos_useLogTables) then
              resultTT(EOS_TABINT_DERIV_DD,g) = &
                   resultTT(EOS_TABINT_DERIV_0,g) * resultTT(EOS_TABINT_DERIV_DD,g) &
                   / speciesDensity
           end if
        end if
     end if
  end do                        ! do g
!
!
!
!


!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabGetSpeciesAnyTableData
