!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_getSpeciesPATableOpacity
!!
!! NAME
!!
!!  op_getSpeciesPATableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesPATableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityPA)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Planck Absorption opacity value for the specified
!!  temperature, density, energy group and species. The interpolation method is
!!  bilinear.
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature
!!  speciesDensity     : The species density
!!  speciesEnergyGroup : The species energy group index
!!  opacityPA          : The value of the determined Planck Absorption opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesPATableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPA           )

  use Driver_interface,  ONLY : Driver_abortFlash

  use op_tabulatedData,  ONLY : op_useLogTables,           &
                                op_nstepsDensityPA,        &
                                op_nstepsTemperaturePA,    &
                                op_tableDensityPA,         &
                                op_tableTemperaturePA,     &
                                op_PlanckAbsorptionTables, &
                                op_species2PATableIndex

  use op_numericsData,   ONLY : one,ten

  implicit none

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityPA

  logical :: lowerBoundary
  logical :: upperBoundary
  logical :: tableBoundary

  integer :: d,t
  integer :: i,j,k,l
  integer :: indexPA
  integer :: nstepsDensity
  integer :: nstepsTemperature

  real :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityPA
  real :: speciesTemperaturePA
  real :: tau,delta
!
!
!   ...Set the handle to the PA tables and associated data arrays.
!
!
  indexPA = op_species2PATableIndex (species)

  if (indexPA == 0) then
      call Driver_abortFlash ('[op_getSpeciesPATableOpacity] ERROR: no handle to PA tables')
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
!        3) the four tabulated PA opacities (o1,o2,o3,o4)
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
  if (op_useLogTables) then
      speciesTemperaturePA = log10 (speciesTemperature)
      speciesDensityPA     = log10 (speciesDensity)
  else
      speciesTemperaturePA = speciesTemperature
      speciesDensityPA     = speciesDensity
  end if

  nstepsTemperature = op_nstepsTemperaturePA (indexPA)

  lowerBoundary = speciesTemperaturePA < op_tableTemperaturePA (1                ,indexPA)
  upperBoundary = speciesTemperaturePA > op_tableTemperaturePA (nstepsTemperature,indexPA)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do t = 2,nstepsTemperature
         if (op_tableTemperaturePA (t,indexPA) >= speciesTemperaturePA) then
             i  = t - 1
             j  = t
             T1 = op_tableTemperaturePA (i,indexPA)
             T2 = op_tableTemperaturePA (j,indexPA)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          i  = 1
          j  = 1
          T1 = speciesTemperaturePA
          T2 = op_tableTemperaturePA (1,indexPA)
      end if

      if (upperBoundary) then
          i  = nstepsTemperature
          j  = nstepsTemperature
          T1 = op_tableTemperaturePA (nstepsTemperature,indexPA)
          T2 = speciesTemperaturePA
      end if
  end if

  nstepsDensity = op_nstepsDensityPA (indexPA)

  lowerBoundary = speciesDensityPA < op_tableDensityPA (1            ,indexPA)
  upperBoundary = speciesDensityPA > op_tableDensityPA (nstepsDensity,indexPA)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do d = 2,nstepsDensity
         if (op_tableDensityPA (d,indexPA) >= speciesDensityPA) then
             k  = d - 1
             l  = d
             D1 = op_tableDensityPA (k,indexPA)
             D2 = op_tableDensityPA (l,indexPA)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          k  = 1
          l  = 1
          D1 = speciesDensityPA
          D2 = op_tableDensityPA (1,indexPA)
      end if

      if (upperBoundary) then
          k  = nstepsDensity
          l  = nstepsDensity
          D1 = op_tableDensityPA (nstepsDensity,indexPA)
          D2 = speciesDensityPA
      end if
  end if

  o1 = op_PlanckAbsorptionTables (i,k,speciesEnergyGroup,indexPA)
  o2 = op_PlanckAbsorptionTables (i,l,speciesEnergyGroup,indexPA)
  o3 = op_PlanckAbsorptionTables (j,k,speciesEnergyGroup,indexPA)
  o4 = op_PlanckAbsorptionTables (j,l,speciesEnergyGroup,indexPA)
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
  tau   = (speciesTemperaturePA - T1) / (T2 - T1)
  delta = (speciesDensityPA     - D1) / (D2 - D1)

  f1 = (one - tau) * (one - delta)
  f2 = delta * (one - tau)
  f3 = tau * (one - delta)
  f4 = delta * tau

  opacityPA = o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4
!
!
!   ...Convert logarithmic form to real form (if needed).
!
!
  if (op_useLogTables) then
      opacityPA = ten ** opacityPA
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesPATableOpacity
