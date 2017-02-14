!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_getSpeciesPETableOpacity
!!
!! NAME
!!
!!  op_getSpeciesPETableOpacity
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesPETableOpacity (integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: speciesEnergyGroup,
!!                                    real    (out) :: opacityPE)
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation a Planck Emission opacity value for the specified
!!  temperature, density, energy group and species. The interpolation method is
!!  bilinear.
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature
!!  speciesDensity     : The species density
!!  speciesEnergyGroup : The species energy group index
!!  opacityPE          : The value of the determined Planck Emission opacity (in cm^2/g)
!!
!!***
subroutine op_getSpeciesPETableOpacity (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityPE           )

  use Driver_interface,  ONLY : Driver_abortFlash

  use op_tabulatedData,  ONLY : op_useLogTables,           &
                                op_nstepsDensityPE,        &
                                op_nstepsTemperaturePE,    &
                                op_tableDensityPE,         &
                                op_tableTemperaturePE,     &
                                op_PlanckEmissionTables,   &
                                op_species2PETableIndex

  use op_numericsData,   ONLY : one,ten

  implicit none
  
  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  real,    intent (out) :: opacityPE

  logical :: lowerBoundary
  logical :: upperBoundary
  logical :: tableBoundary

  integer :: d,t
  integer :: i,j,k,l
  integer :: indexPE
  integer :: nstepsDensity
  integer :: nstepsTemperature

  real :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityPE
  real :: speciesTemperaturePE
  real :: tau,delta
!
!
!   ...Set the handle to the PE tables and associated data arrays.
!
!
  indexPE = op_species2PETableIndex (species)

  if (indexPE == 0) then
      call Driver_abortFlash ('[op_getSpeciesPETableOpacity] ERROR: no handle to PE tables')
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
!        3) the four tabulated PE opacities (o1,o2,o3,o4)
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
      speciesTemperaturePE = log10 (speciesTemperature)
      speciesDensityPE     = log10 (speciesDensity)
  else
      speciesTemperaturePE = speciesTemperature
      speciesDensityPE     = speciesDensity
  end if

  nstepsTemperature = op_nstepsTemperaturePE (indexPE)

  lowerBoundary = speciesTemperaturePE < op_tableTemperaturePE (1                ,indexPE)
  upperBoundary = speciesTemperaturePE > op_tableTemperaturePE (nstepsTemperature,indexPE)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do t = 2,nstepsTemperature
         if (op_tableTemperaturePE (t,indexPE) >= speciesTemperaturePE) then
             i  = t - 1
             j  = t
             T1 = op_tableTemperaturePE (i,indexPE)
             T2 = op_tableTemperaturePE (j,indexPE)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          i  = 1
          j  = 1
          T1 = speciesTemperaturePE
          T2 = op_tableTemperaturePE (1,indexPE)
      end if

      if (upperBoundary) then
          i  = nstepsTemperature
          j  = nstepsTemperature
          T1 = op_tableTemperaturePE (nstepsTemperature,indexPE)
          T2 = speciesTemperaturePE
      end if
  end if

  nstepsDensity = op_nstepsDensityPE (indexPE)

  lowerBoundary = speciesDensityPE < op_tableDensityPE (1            ,indexPE)
  upperBoundary = speciesDensityPE > op_tableDensityPE (nstepsDensity,indexPE)
  tableBoundary = (.not.lowerBoundary) .and. (.not.upperBoundary)

  if (tableBoundary) then
      do d = 2,nstepsDensity
         if (op_tableDensityPE (d,indexPE) >= speciesDensityPE) then
             k  = d - 1
             l  = d
             D1 = op_tableDensityPE (k,indexPE)
             D2 = op_tableDensityPE (l,indexPE)
             exit
         end if
      end do
  else
      if (lowerBoundary) then
          k  = 1
          l  = 1
          D1 = speciesDensityPE
          D2 = op_tableDensityPE (1,indexPE)
      end if

      if (upperBoundary) then
          k  = nstepsDensity
          l  = nstepsDensity
          D1 = op_tableDensityPE (nstepsDensity,indexPE)
          D2 = speciesDensityPE
      end if
  end if

  o1 = op_PlanckEmissionTables (i,k,speciesEnergyGroup,indexPE)
  o2 = op_PlanckEmissionTables (i,l,speciesEnergyGroup,indexPE)
  o3 = op_PlanckEmissionTables (j,k,speciesEnergyGroup,indexPE)
  o4 = op_PlanckEmissionTables (j,l,speciesEnergyGroup,indexPE)
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
  tau   = (speciesTemperaturePE - T1) / (T2 - T1)
  delta = (speciesDensityPE     - D1) / (D2 - D1)

  f1 = (one - tau) * (one - delta)
  f2 = delta * (one - tau)
  f3 = tau * (one - delta)
  f4 = delta * tau

  opacityPE = o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4
!
!
!   ...Convert logarithmic form to real form (if needed).
!
!
  if (op_useLogTables) then
      opacityPE = ten ** opacityPE
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesPETableOpacity
