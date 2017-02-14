!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_getOpalTableOpacity
!!
!! NAME
!!
!!  op_getOpalTableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getOpalTableOpacity (integer (in)  :: tempRange,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityRO)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Rosseland opacity value for
!!  the specified energy group and tempRange. The interpolation method
!!  is bilinear.
!!
!! ARGUMENTS
!!
!!  tempRange            : The tempRange index
!!  speciesTemperature : The species temperature
!!  speciesDensity     : The species density
!!  speciesEnergyGroup : The species energy group index
!!  opacityRO          : The value of the determined Rosseland opacity (in cm^2/g)
!!
!!***
subroutine op_getOpalTableOpacity (mfH, tempRange,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityRO           )

  use Driver_interface,  ONLY : Driver_abortFlash

  use op_opalData,  ONLY : op_useLogTables,           &
                                opT_varTableGroupPT,       &
                                op_opalAllTab

  implicit none

  integer, intent (in)  :: mfH
  integer, intent (in)  :: tempRange
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityRO

  real, parameter :: one   = 1.0
  real, parameter :: ten   = 10.0

  type(opT_varTableGroupPT),pointer :: thisTypeTable
  real,    pointer :: thisTypeDensity (:)
  real,    pointer :: thisTypeTemperature (:)

  logical :: lowerBoundary
  logical :: upperBoundary
  logical :: withinBoundary

  integer :: d,t
  integer :: i,j,k,l
  integer :: nstepsDensity
  integer :: nstepsTemperature

  real :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityRO
  real :: speciesTemperatureRO
  real :: logR, R
  real :: tau,delta
!
!
!   ...Set the handle to the RO tables and associated data arrays.
!
!

  thisTypeTable => op_opalAllTab(mfH,tempRange)%tg

  if (associated(thisTypeTable)) then
#ifdef DEBUG_OPACITY
     print*,'op_getOpalTableOpacity: mfH,temprange are', mfH,tempRange,' ,ok'
#endif
  else
     print*,'op_getOpalTableOpacity: mfH,temprange are', mfH,tempRange
     call Driver_abortFlash ('[op_getOpalTableOpacity] ERROR: thisTypeTable not associated; tg NULL?')
  end if

  if (associated(thisTypeTable%table)) then
#ifdef DEBUG_OPACITY
     print*,'op_getOpalTableOpacity: thisTypeTable%table ok'
#endif
  else
     print*,'op_getOpalTableOpacity: mfH,temprange are', mfH,tempRange
     call Driver_abortFlash ('[op_getOpalTableOpacity] ERROR: thisTypeTable%table not associated; tg%table NULL?')
  end if

  if (associated(thisTypeTable%table%table)) then
#ifdef DEBUG_OPACITY
     print*,'op_getOpalTableOpacity: thisTypeTable%table%table ok'
#endif
  else
     print*,'op_getOpalTableOpacity: mfH,temprange are', mfH,tempRange
     call Driver_abortFlash ('[op_getOpalTableOpacity] ERROR: table%table not associated; tg%table%table NULL?')
  end if

  thisTypeDensity => thisTypeTable%td%Densities
  thisTypeTemperature => thisTypeTable%td%Temperatures

!
!
!   ...Get the current temperature and density and find:
!
!        1) the temperature/density quadrant (T1,T2,D1,D2) boundary containing
!           containing the temperature and density (x)
!
!        2) the table index quadrant (i,j,k,l) of the boundary, as indicated
!           in the figure below.
!
!        3) the four tabulated RO opacities (o1,o2,o3,o4)
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
!      In case the temperature and/or density lay outside
!      the tabulated boundaries, take the corresponding boundary values
!      of the tables. The criteria as to when the temperature and
!      density belong within the boundary are:
!
!                         T1 =< speciesTemperature =< T2
!                         D1 =<   speciesDensity   =< D2
!
!
  if (op_useLogTables) then
      speciesTemperatureRO = log10 (speciesTemperature)
      speciesDensityRO     = log10 (speciesDensity)
!!$      logR = logRho - 3*logT + 18
      logR = speciesDensityRO - 3.0*speciesTemperatureRO + 18.0
      R = logR
  else
      speciesTemperatureRO = speciesTemperature
      speciesDensityRO     = speciesDensity
!!$      R = Rho / T**3 * 1.E18
      R = speciesDensity / speciesTemperature**3 * 1.E18
  end if

  nstepsTemperature = thisTypeTable%td%ntemp

  lowerBoundary = speciesTemperatureRO < thisTypeTemperature (1            )
  if (lowerBoundary) then
     upperBoundary = .FALSE.
  else
     upperBoundary = speciesTemperatureRO > thisTypeTemperature (nstepsTemperature) .OR. nstepsTemperature == 1
  end if
  withinBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (withinBoundary) then
      do t = 2,nstepsTemperature
         if (thisTypeTemperature (t) >= speciesTemperatureRO) then
             i  = t - 1
             j  = t
             T1 = thisTypeTemperature (i)
             T2 = thisTypeTemperature (j)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          i  = 1
          j  = 1
          T1 = speciesTemperatureRO
          T2 = thisTypeTemperature (1)
      end if

      if (upperBoundary) then
          i  = nstepsTemperature
          j  = nstepsTemperature
          T1 = thisTypeTemperature (nstepsTemperature)
          T2 = speciesTemperatureRO
      end if
  end if

  nstepsDensity = thisTypeTable%td%ndens

  lowerBoundary =        R         < thisTypeDensity (1            )
  if (lowerBoundary) then
     upperBoundary = .FALSE.
  else
     upperBoundary =        R         > thisTypeDensity (nstepsDensity) .OR. nstepsDensity == 1
  end if
  withinBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)
#ifdef DEBUG_OPACITY
  print*,'speciesDensity{,RO}=',speciesDensity,speciesDensityRO,&
       thisTypeDensity (1            ),thisTypeDensity (nstepsDensity)
  print*,'nstepsDensity=',nstepsDensity,lowerBoundary,upperBoundary,withinBoundary
#endif
  if (withinBoundary) then
      do d = 2,nstepsDensity
         if (thisTypeDensity (d) >=        R        ) then
             k  = d - 1
             l  = d
             D1 = thisTypeDensity (k)
             D2 = thisTypeDensity (l)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          k  = 1
          l  = 1
          D1 =        R        
          D2 = thisTypeDensity (1)
      end if

      if (upperBoundary) then
          k  = nstepsDensity
          l  = nstepsDensity
          D1 = thisTypeDensity (nstepsDensity)
          D2 =        R        
      end if
  end if

#ifdef DEBUG_OPACITY
  print*,'SHAPE(thisTypeTable%table%table)=',SHAPE(thisTypeTable%table%table)
  print*,'[op_getOpalTableOpacity] i,j,k,l=',i,j,k,l
#endif
  o1 = thisTypeTable%table%table(i,k)
  o2 = thisTypeTable%table%table(i,l)
  o3 = thisTypeTable%table%table(j,k)
  o4 = thisTypeTable%table%table(j,l)
!
!
!   ...Do the bilinear interpolation:
!
!               opacity =   o1 * [(1-tau)*(1-delta)]
!                         + o2 * [delta*(1-tau)]
!                         + o3 * [tau*(1-delta)]
!                         + o4 * [delta*tau]
!
!
  tau   = (speciesTemperatureRO - T1) / (T2 - T1)
  delta = (       R             - D1) / (D2 - D1)

  f1 = (one - tau) * (one - delta)
  f2 = delta * (one - tau)
  f3 = tau * (one - delta)
  f4 = delta * tau

  opacityRO = o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4
!
!
!   ...Convert logarithmic form to real form (if needed).
!
!
  if (op_useLogTables) then
      opacityRO = ten ** opacityRO
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getOpalTableOpacity
