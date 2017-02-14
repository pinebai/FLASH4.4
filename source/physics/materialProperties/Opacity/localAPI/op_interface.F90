!!****ih* source/physics/materialProperties/Opacity/localAPI/op_interface
!!
!! NAME
!!
!!  op_interface
!!
!! SYNOPSIS
!!
!!  op_interface ()
!!
!! DESCRIPTION
!!
!!  This is the interface file for the local API in the Opacity Unit.
!!
!!***

module op_interface
  
  interface
     subroutine op_BiggsGroupOpacity (opacityKind,indexElement,Temperature,Elower,Eupper,Opacity)
       character (len=9), intent (in)  :: opacityKind
       integer,           intent (in)  :: indexElement
       real,              intent (in)  :: Temperature
       real,              intent (in)  :: Elower
       real,              intent (in)  :: Eupper
       real,              intent (out) :: Opacity
     end subroutine op_BiggsGroupOpacity
  end interface

  interface
     subroutine op_BiggsPlanckGroupIntegrate (nSec,                                &
                                              A1Sec,A2Sec,A3Sec,A4Sec,             &
                                              intLimits,                           &
                                                            rescaleBase10Exponent, &
                                                            BiggsPlanckIntegral,   &
                                                            PlanckIntegral         )
       integer, intent (in)  :: nSec
       real,    intent (in)  :: A1Sec     (1:nSec)
       real,    intent (in)  :: A2Sec     (1:nSec)
       real,    intent (in)  :: A3Sec     (1:nSec)
       real,    intent (in)  :: A4Sec     (1:nSec)
       real,    intent (in)  :: intLimits (1:nSec+1)
       real,    intent (out) :: rescaleBase10Exponent
       real,    intent (out) :: BiggsPlanckIntegral
       real,    intent (out) :: PlanckIntegral
     end subroutine op_BiggsPlanckGroupIntegrate
  end interface

  interface
     subroutine op_BiggsRosslndGroupIntegrate (nSec,                                &
                                               A1Sec,A2Sec,A3Sec,A4Sec,             &
                                               intLimits,                           &
                                                             rescaleBase10Exponent, &
                                                             BiggsRosslndIntegral,  &
                                                             RosslndIntegral        )
       integer, intent (in)  :: nSec
       real,    intent (in)  :: A1Sec     (1:nSec)
       real,    intent (in)  :: A2Sec     (1:nSec)
       real,    intent (in)  :: A3Sec     (1:nSec)
       real,    intent (in)  :: A4Sec     (1:nSec)
       real,    intent (in)  :: intLimits (1:nSec+1)
       real,    intent (out) :: rescaleBase10Exponent
       real,    intent (out) :: BiggsRosslndIntegral
       real,    intent (out) :: RosslndIntegral
     end subroutine op_BiggsRosslndGroupIntegrate
  end interface

  interface
     subroutine op_browseIonmixTables (tableName,                        &
                                       needPATable,                      &
                                       needPETable,                      &
                                       needROTable,                      &
                                                    nstepsDensityPA,     &
                                                    nstepsDensityPE,     &
                                                    nstepsDensityRO,     &
                                                    nstepsTemperaturePA, &
                                                    nstepsTemperaturePE, &
                                                    nstepsTemperatureRO  )

       character (len=80), intent (in)  :: tableName
       logical,            intent (in)  :: needPATable
       logical,            intent (in)  :: needPETable
       logical,            intent (in)  :: needROTable
       integer,            intent (out) :: nstepsDensityPA
       integer,            intent (out) :: nstepsDensityPE
       integer,            intent (out) :: nstepsDensityRO
       integer,            intent (out) :: nstepsTemperaturePA
       integer,            intent (out) :: nstepsTemperaturePE
       integer,            intent (out) :: nstepsTemperatureRO
     end subroutine op_browseIonmixTables
  end interface

  interface
     subroutine op_browseTables (tableKind,                        &
                                 tableName,                        &
                                 needPATable,                      &
                                 needPETable,                      &
                                 needROTable,                      &
                                              nstepsDensityPA,     &
                                              nstepsDensityPE,     &
                                              nstepsDensityRO,     &
                                              nstepsTemperaturePA, &
                                              nstepsTemperaturePE, &
                                              nstepsTemperatureRO  )

       character (len=80), intent (in)  :: tableKind
       character (len=80), intent (in)  :: tableName
       logical,            intent (in)  :: needPATable
       logical,            intent (in)  :: needPETable
       logical,            intent (in)  :: needROTable
       integer,            intent (out) :: nstepsDensityPA
       integer,            intent (out) :: nstepsDensityPE
       integer,            intent (out) :: nstepsDensityRO
       integer,            intent (out) :: nstepsTemperaturePA
       integer,            intent (out) :: nstepsTemperaturePE
       integer,            intent (out) :: nstepsTemperatureRO
     end subroutine op_browseTables
  end interface

  interface
     subroutine op_browseOpalTable (tableName,                        &
                                                    nstepsDensity,    &
                                                    nstepsTemperature  )

       character (len=80), intent (in)  :: tableName
       integer,            intent (out) :: nstepsDensity
       integer,            intent (out) :: nstepsTemperature
     end subroutine op_browseOpalTable
  end interface

  interface
     subroutine op_calculateLowTempOpacities ()
     end subroutine op_calculateLowTempOpacities
  end interface

  interface
     subroutine op_computeIonNumberDensities ()
     end subroutine op_computeIonNumberDensities
  end interface

  interface
     subroutine op_computeMultispeciesOpacities (opacityAbsorption, &
                                                 opacityEmission,   &
                                                 opacityTransport   )

       real,    intent (out) :: opacityAbsorption
       real,    intent (out) :: opacityEmission
       real,    intent (out) :: opacityTransport
     end subroutine op_computeMultispeciesOpacities
  end interface

  interface
     subroutine op_finalizeConstant ()
     end subroutine op_finalizeConstant
  end interface

  interface
     subroutine op_finalizeConstcm2g ()
     end subroutine op_finalizeConstcm2g
  end interface

  interface
     subroutine op_finalizeIntegrate ()
     end subroutine op_finalizeIntegrate
  end interface

  interface
     subroutine op_finalizeLowTemp ()
     end subroutine op_finalizeLowTemp
  end interface

  interface
     subroutine op_finalizeNumerics ()
     end subroutine op_finalizeNumerics
  end interface

  interface
     subroutine op_finalizeTabulated ()
     end subroutine op_finalizeTabulated
  end interface

  interface
     subroutine op_generateRootsWeights (n,nV,MomZero,V,Diag,Offd,Roots,Weights)
       integer, intent (in)    :: n
       integer, intent (in)    :: nV
       real,    intent (in)    :: MomZero
       real,    intent (inout) :: V       (1:nV)
       real,    intent (inout) :: Diag    (1:n)
       real,    intent (inout) :: Offd    (1:n)
       real,    intent (out)   :: Roots   (1:n)
       real,    intent (out)   :: Weights (1:n)
     end subroutine op_generateRootsWeights
  end interface

  interface
     subroutine op_getSpeciesConstantOpacities (species, dens)
       integer, intent (in) :: species
       real,    intent (in) :: dens
     end subroutine op_getSpeciesConstantOpacities
  end interface

  interface
     subroutine op_getSpeciesConstcm2gOpacities (species)
       integer, intent (in) :: species
     end subroutine op_getSpeciesConstcm2gOpacities
  end interface

  interface
     subroutine op_getSpeciesLowTempOpacities (species,            &
                                               speciesTemperature, &
                                               speciesEnergyGroup  )
       integer, intent (in)  :: species
       integer, intent (in)  :: speciesEnergyGroup
       real,    intent (in)  :: speciesTemperature
     end subroutine op_getSpeciesLowTempOpacities
  end interface

  interface
     subroutine op_getSpeciesPATableOpacity (species,            &
                                             speciesTemperature, &
                                             speciesDensity,     &
                                             speciesEnergyGroup, &
                                             opacityPA           )
       integer, intent (in)  :: species
       integer, intent (in)  :: speciesEnergyGroup
       real,    intent (in)  :: speciesTemperature
       real,    intent (in)  :: speciesDensity
       real,    intent (out) :: opacityPA
     end subroutine op_getSpeciesPATableOpacity
  end interface

  interface
     subroutine op_getSpeciesPETableOpacity (species,            &
                                             speciesTemperature, &
                                             speciesDensity,     &
                                             speciesEnergyGroup, &
                                             opacityPE           )
       integer, intent (in)  :: species
       integer, intent (in)  :: speciesEnergyGroup
       real,    intent (in)  :: speciesTemperature
       real,    intent (in)  :: speciesDensity
       real,    intent (out) :: opacityPE
     end subroutine op_getSpeciesPETableOpacity
  end interface

  interface
     subroutine op_getSpeciesROTableOpacity (species,            &
                                             speciesTemperature, &
                                             speciesDensity,     &
                                             speciesEnergyGroup, &
                                             opacityRO           )
       integer, intent (in)  :: species
       integer, intent (in)  :: speciesEnergyGroup
       real,    intent (in)  :: speciesTemperature
       real,    intent (in)  :: speciesDensity
       real,    intent (out) :: opacityRO
     end subroutine op_getSpeciesROTableOpacity
  end interface

  interface
     subroutine op_getSpeciesTabulatedOpacities (species,            &
                                                 speciesTemperature, &
                                                 speciesDensity,     &
                                                 speciesEnergyGroup, &
                                                 needPATable,        &
                                                 needPETable,        &
                                                 needROTable         )
       logical, intent (in) :: needPATable
       logical, intent (in) :: needPETable
       logical, intent (in) :: needROTable
       integer, intent (in) :: species
       integer, intent (in) :: speciesEnergyGroup
       real,    intent (in) :: speciesTemperature
       real,    intent (in) :: speciesDensity
     end subroutine op_getSpeciesTabulatedOpacities
  end interface

  interface
     subroutine op_getOpalTableOpacity (mfH, tempRange,     &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        speciesEnergyGroup, &
                                        opacityRO           )
       implicit none
       integer, intent (in)  :: mfH, tempRange
       integer, intent (in)  :: speciesEnergyGroup
       real,    intent (in)  :: speciesTemperature
       real,    intent (in)  :: speciesDensity
       real,    intent (out) :: opacityRO
     end subroutine op_getOpalTableOpacity
  end interface

  interface
     subroutine op_initConstant ()
     end subroutine op_initConstant
  end interface

  interface
     subroutine op_initConstcm2g ()
     end subroutine op_initConstcm2g
  end interface

  interface
     subroutine op_initIntegrate ()
     end subroutine op_initIntegrate
  end interface

  interface
     subroutine op_initLowTemp ()
     end subroutine op_initLowTemp
  end interface

  interface
     subroutine op_initNumerics ()
     end subroutine op_initNumerics
  end interface

  interface
     subroutine op_initTabulated ()
     end subroutine op_initTabulated
  end interface

  interface
     subroutine op_KleinGroupOpacity (opacityKind,indexElement,Temperature,Elower,Eupper,Opacity)
       character (len=9), intent (in)  :: opacityKind
       integer,           intent (in)  :: indexElement
       real,              intent (in)  :: Temperature
       real,              intent (in)  :: Elower
       real,              intent (in)  :: Eupper
       real,              intent (out) :: Opacity
     end subroutine op_KleinGroupOpacity
  end interface

  interface
     subroutine op_KleinPlanckGroupIntegrate (KNintegralFactor,              &
                                              B,                             &
                                              R,S,                           &
                                                      rescaleBase10Exponent, &
                                                      KleinPlanckIntegral,   &
                                                      PlanckIntegral         )
       real,    intent (in)  :: KNintegralFactor
       real,    intent (in)  :: B
       real,    intent (in)  :: R,S
       real,    intent (out) :: rescaleBase10Exponent
       real,    intent (out) :: KleinPlanckIntegral
       real,    intent (out) :: PlanckIntegral
     end subroutine op_KleinPlanckGroupIntegrate
  end interface

  interface
     subroutine op_KleinRosslndGroupIntegrate (KNintegralFactor,              &
                                               B,                             &
                                               R,S,                           &
                                                       rescaleBase10Exponent, &
                                                       KleinRosslndIntegral,  &
                                                       RosslndIntegral        )
       real,    intent (in)  :: KNintegralFactor
       real,    intent (in)  :: B
       real,    intent (in)  :: R,S
       real,    intent (out) :: rescaleBase10Exponent
       real,    intent (out) :: KleinRosslndIntegral
       real,    intent (out) :: RosslndIntegral
     end subroutine op_KleinRosslndGroupIntegrate
  end interface

  interface
     subroutine op_LaguerreCoefficients (n,alpha,LagA,LagB)
       integer, intent (in)  :: n
       real,    intent (in)  :: alpha
       real,    intent (out) :: LagA (1:n)
       real,    intent (out) :: LagB (1:n)
     end subroutine op_LaguerreCoefficients
  end interface

  interface
     subroutine op_LaguerreMoments (n,nLagT,alpha,beta,T,LagT,MomZero,Mom)
       integer, intent (in)    :: n
       integer, intent (in)    :: nLagT
       real,    intent (in)    :: alpha
       real,    intent (in)    :: beta
       real,    intent (in)    :: T
       real,    intent (out)   :: MomZero
       real,    intent (inout) :: LagT (1:nLagT)
       real,    intent (out)   :: Mom  (1:n)
     end subroutine op_LaguerreMoments
  end interface

  interface
     subroutine op_LaguerrePolynomials (n,nLag,alpha,X,Lag)
       integer, intent (in)  :: n
       integer, intent (in)  :: nLag
       real,    intent (in)  :: alpha
       real,    intent (in)  :: X
       real,    intent (out) :: Lag (1:nLag)
     end subroutine op_LaguerrePolynomials
  end interface

  interface
     subroutine op_LaguerreQuadratureRule (nRoots,beta,T,Roots,Weights)
       integer, intent (in)  :: nRoots
       real,    intent (in)  :: beta
       real,    intent (in)  :: T
       real,    intent (out) :: Roots   (1:nRoots)
       real,    intent (out) :: Weights (1:nRoots)
     end subroutine op_LaguerreQuadratureRule
  end interface

  interface
     subroutine op_LaguerreZeroMoment (beta,R,S,MomZero)
       real, intent (in)  :: beta
       real, intent (in)  :: R,S
       real, intent (out) :: MomZero
     end subroutine op_LaguerreZeroMoment
  end interface

  interface
     subroutine op_OrthoPolyCoefficients (n,m,nRow,Mom,A,B,Row1,Row2,OrthoA,OrthoB)
       integer, intent (in)    :: n
       integer, intent (in)    :: m
       integer, intent (in)    :: nRow
       real,    intent (in)    :: Mom    (1:m)
       real,    intent (in)    :: A      (1:m)
       real,    intent (in)    :: B      (1:m)
       real,    intent (inout) :: Row1   (1:nRow)
       real,    intent (inout) :: Row2   (1:nRow)
       real,    intent (out)   :: OrthoA (1:n)
       real,    intent (out)   :: OrthoB (1:n)
     end subroutine op_OrthoPolyCoefficients
  end interface

  interface
     subroutine op_readBiggs1971xraydatFile ()
     end subroutine op_readBiggs1971xraydatFile
  end interface

  interface
     subroutine op_readIonmixTables (tableName,   &
                                     needPATable, &
                                     needPETable, &
                                     needROTable, &
                                     indexPA,     &
                                     indexPE,     &
                                     indexRO      )

       character (len=80), intent (in) :: tableName
       logical,            intent (in) :: needPATable
       logical,            intent (in) :: needPETable
       logical,            intent (in) :: needROTable
       integer,            intent (in) :: indexPA
       integer,            intent (in) :: indexPE
       integer,            intent (in) :: indexRO
     end subroutine op_readIonmixTables
  end interface

  interface
     subroutine op_readIonmix4Tables (tableName,   &
                                      needPATable, &
                                      needPETable, &
                                      needROTable, &
                                      indexPA,     &
                                      indexPE,     &
                                      indexRO      )

       character (len=80), intent (in) :: tableName
       logical,            intent (in) :: needPATable
       logical,            intent (in) :: needPETable
       logical,            intent (in) :: needROTable
       integer,            intent (in) :: indexPA
       integer,            intent (in) :: indexPE
       integer,            intent (in) :: indexRO
     end subroutine op_readIonmix4Tables
  end interface


  interface
     subroutine op_readTables (tableKind,   &
                               tableName,   &
                               needPATable, &
                               needPETable, &
                               needROTable, &
                               indexPA,     &
                               indexPE,     &
                               indexRO      )

       character (len=80), intent (in) :: tableKind
       character (len=80), intent (in) :: tableName
       logical,            intent (in) :: needPATable
       logical,            intent (in) :: needPETable
       logical,            intent (in) :: needROTable
       integer,            intent (in) :: indexPA
       integer,            intent (in) :: indexPE
       integer,            intent (in) :: indexRO
     end subroutine op_readTables
  end interface


  interface
     subroutine op_setAtomNames ()
     end subroutine op_setAtomNames
  end interface

  interface
     subroutine op_setAtomWeights ()
     end subroutine op_setAtomWeights
  end interface

  interface
     subroutine op_setEnergyGroupBoundaries ()
     end subroutine op_setEnergyGroupBoundaries
  end interface

  interface
     subroutine op_setLowTemperatureBoundaries ()
     end subroutine op_setLowTemperatureBoundaries
  end interface

  interface
     subroutine op_setPEarrayJmax ()
     end subroutine op_setPEarrayJmax
  end interface

  interface
     subroutine op_setPEcoeffsAij4 ()
     end subroutine op_setPEcoeffsAij4
  end interface

  interface
     subroutine op_setPEcoeffsAij4atomsZ100 ()
     end subroutine op_setPEcoeffsAij4atomsZ100
  end interface

  interface
     subroutine op_setPEcoeffsAij4atomsZ43 ()
     end subroutine op_setPEcoeffsAij4atomsZ43
  end interface

  interface
     subroutine op_setPEcoeffsAij4atomsZ73 ()
     end subroutine op_setPEcoeffsAij4atomsZ73
  end interface

  interface
     subroutine op_setPEenergyRange ()
     end subroutine op_setPEenergyRange
  end interface

  interface
     subroutine op_setPEenergyRangeAtomsZ100 ()
     end subroutine op_setPEenergyRangeAtomsZ100
  end interface

  interface
     subroutine op_setPEenergyRangeAtomsZ43 ()
     end subroutine op_setPEenergyRangeAtomsZ43
  end interface

  interface
     subroutine op_setPEenergyRangeAtomsZ73 ()
     end subroutine op_setPEenergyRangeAtomsZ73
  end interface

  interface
     subroutine op_setSpeciesElementsData ()
     end subroutine op_setSpeciesElementsData
  end interface

  interface
     subroutine op_shJacobi3TermMoments (n,maxM,p,beta,T,X,Y,Sinv,MomZero,Mom)
       integer, intent (in)    :: n
       integer, intent (in)    :: maxM
       real,    intent (in)    :: p
       real,    intent (in)    :: beta
       real,    intent (in)    :: T
       real,    intent (out)   :: MomZero
       real,    intent (inout) :: X    (1:maxM)
       real,    intent (inout) :: Y    (1:maxM)
       real,    intent (inout) :: Sinv (1:maxM)
       real,    intent (out)   :: Mom  (1:n)
     end subroutine op_shJacobi3TermMoments
  end interface

  interface
     subroutine op_shJacobiCoefficients (n,p,q,JacA,JacB)
       integer, intent (in)  :: n
       real,    intent (in)  :: p,q
       real,    intent (out) :: JacA (1:n)
       real,    intent (out) :: JacB (1:n)
     end subroutine op_shJacobiCoefficients
  end interface

  interface
     subroutine op_writeAtomPEopacity2file (Z,Elower,Eupper,nPoints)
       integer, intent (in) :: nPoints
       integer, intent (in) :: Z
       real,    intent (in) :: Elower
       real,    intent (in) :: Eupper
     end subroutine op_writeAtomPEopacity2file
  end interface

  interface
     subroutine op_writeCompletePEdata ()
     end subroutine op_writeCompletePEdata
  end interface

  interface
     subroutine op_writeConstants ()
     end subroutine op_writeConstants
  end interface

  interface
     subroutine op_writeConstcm2g ()
     end subroutine op_writeConstcm2g
  end interface

  interface
     subroutine op_writeElementsPEdata ()
     end subroutine op_writeElementsPEdata
  end interface

  interface
     subroutine op_writeLowTempTables ()
     end subroutine op_writeLowTempTables
  end interface

  interface
     subroutine op_writeNumerics ()
     end subroutine op_writeNumerics
  end interface

  interface
     subroutine op_writeOpacity ()
     end subroutine op_writeOpacity
  end interface

  interface
     subroutine op_writeQuadratureData (nRoots,nMoments,              &
                                        Ttiny,Tsmall,Tlarge,Thuge,T,  &
                                        p,q,beta,MomZero,Roots,Weights)
       integer, intent (in) :: nRoots
       integer, intent (in) :: nMoments
       logical, intent (in) :: Ttiny
       logical, intent (in) :: Tsmall
       logical, intent (in) :: Tlarge
       logical, intent (in) :: Thuge
       real,    intent (in) :: T
       real,    intent (in) :: p,q
       real,    intent (in) :: beta
       real,    intent (in) :: MomZero
       real,    intent (in) :: Roots   (1:nRoots)
       real,    intent (in) :: Weights (1:nRoots)
     end subroutine op_writeQuadratureData
  end interface

  interface
     subroutine op_writeSpeciesPATable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine op_writeSpeciesPATable
  end interface

  interface
     subroutine op_writeSpeciesPETable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine op_writeSpeciesPETable
  end interface

  interface
     subroutine op_writeSpeciesROTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine op_writeSpeciesROTable
  end interface

  interface
     subroutine op_writeSpeciesTables (fileUnit,    &
                                       species,     &
                                       needPATable, &
                                       needPETable, &
                                       needROTable  )

       logical, intent (in) :: needPATable
       logical, intent (in) :: needPETable
       logical, intent (in) :: needROTable
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine op_writeSpeciesTables
  end interface

  interface
     subroutine op_writeTables ()
     end subroutine op_writeTables
  end interface

end module op_interface
