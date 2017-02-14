!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_getSpeciesROTableOpacity
!!
!! NAME
!!
!!  op_getSpeciesROTableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesROTableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityRO)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Rosseland opacity value for
!!  the specified energy group and species. The interpolation method
!!  is bilinear.
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature
!!  speciesDensity     : The species density
!!  speciesEnergyGroup : The species energy group index
!!  opacityRO          : The value of the determined Rosseland opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesROTableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityRO           )

  use Driver_interface,  ONLY : Driver_abortFlash

  use op_tabulatedData,  ONLY : op_useLogTables,           &
                                op_nstepsDensityRO,        &
                                op_nstepsTemperatureRO,    &
                                op_tableDensityRO,         &
                                op_tableTemperatureRO,     &
                                op_RosselandTables,        &
                                op_species2ROTableIndex

  use op_numericsData,   ONLY : one,ten

  implicit none

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityRO

  logical :: lowerBoundary
  logical :: upperBoundary
  logical :: tableBoundary

  integer :: d,t
  integer :: i,j,k,l
  integer :: indexRO
  integer :: nstepsDensity
  integer :: nstepsTemperature

  real :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityRO
  real :: speciesTemperatureRO
  real :: tau,delta
!
!
!   ...Set the handle to the RO tables and associated data arrays.
!
!
  indexRO = op_species2ROTableIndex (species)

  if (indexRO == 0) then
      call Driver_abortFlash ('[op_getSpeciesROTableOpacity] ERROR: no handle to RO tables')
  end if
!
!
!   ...Get the current temperature and density of the species and find:
!
!        1) the temperature/density quadrant (T1,T2,D1,D2) boundary containing
!           containing the species's temperature and density (x)
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
      speciesTemperatureRO = log10 (speciesTemperature)
      speciesDensityRO     = log10 (speciesDensity)
  else
      speciesTemperatureRO = speciesTemperature
      speciesDensityRO     = speciesDensity
  end if

  nstepsTemperature = op_nstepsTemperatureRO (indexRO)

  lowerBoundary = speciesTemperatureRO < op_tableTemperatureRO (1                ,indexRO)
  upperBoundary = speciesTemperatureRO > op_tableTemperatureRO (nstepsTemperature,indexRO)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do t = 2,nstepsTemperature
         if (op_tableTemperatureRO (t,indexRO) >= speciesTemperatureRO) then
             i  = t - 1
             j  = t
             T1 = op_tableTemperatureRO (i,indexRO)
             T2 = op_tableTemperatureRO (j,indexRO)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          i  = 1
          j  = 1
          T1 = speciesTemperatureRO
          T2 = op_tableTemperatureRO (1,indexRO)
      end if

      if (upperBoundary) then
          i  = nstepsTemperature
          j  = nstepsTemperature
          T1 = op_tableTemperatureRO (nstepsTemperature,indexRO)
          T2 = speciesTemperatureRO
      end if
  end if

  nstepsDensity = op_nstepsDensityRO (indexRO)

  lowerBoundary = speciesDensityRO < op_tableDensityRO (1            ,indexRO)
  upperBoundary = speciesDensityRO > op_tableDensityRO (nstepsDensity,indexRO)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do d = 2,nstepsDensity
         if (op_tableDensityRO (d,indexRO) >= speciesDensityRO) then
             k  = d - 1
             l  = d
             D1 = op_tableDensityRO (k,indexRO)
             D2 = op_tableDensityRO (l,indexRO)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          k  = 1
          l  = 1
          D1 = speciesDensityRO
          D2 = op_tableDensityRO (1,indexRO)
      end if

      if (upperBoundary) then
          k  = nstepsDensity
          l  = nstepsDensity
          D1 = op_tableDensityRO (nstepsDensity,indexRO)
          D2 = speciesDensityRO
      end if
  end if

  o1 = op_RosselandTables (i,k,speciesEnergyGroup,indexRO)
  o2 = op_RosselandTables (i,l,speciesEnergyGroup,indexRO)
  o3 = op_RosselandTables (j,k,speciesEnergyGroup,indexRO)
  o4 = op_RosselandTables (j,l,speciesEnergyGroup,indexRO)
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
  delta = (speciesDensityRO     - D1) / (D2 - D1)

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
end subroutine op_getSpeciesROTableOpacity
