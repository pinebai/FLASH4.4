!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/EnergyDeposition_data
!!
!! NAME
!!
!!  EnergyDeposition_data
!!
!! SYNOPSIS
!!
!!  use EnergyDeposition_data
!!  
!! DESCRIPTION
!!
!!  Data module for EnergyDeposition
!!  --------------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (C) means data that is (C)alculated internally by the laser code
!!
!!           Nele = number of electrons
!!           Tele = electron temperature
!!
!!   ed_Avogadro                    (G) : Avogadro's constant (local copy)
!!   ed_badTorsionAxis              (P) : An internal criterion for having chosen the wrong torsion axis
!!   ed_baseName                    (R) : The simulation base name
!!   ed_beams                     (R,C) : An array of special type, containing all the beams info (see below)
!!   ed_beamsAreSetup               (C) : Logical indicator showing if the beams have been set up
!!   ed_Boltzmann                   (G) : Boltzmann's constant (local copy)
!!   ed_cellCenters                 (C) : Holds the center coordinates for all cells within a block
!!   ed_cellCubicNele               (C) : Holds the Nele cubic expansion coefficients for all cells within a block
!!   ed_cellCubicTele               (C) : Holds the Tele cubic expansion coefficients for all cells within a block
!!   ed_cellDensity                 (C) : Holds the mass density for all cells within a block
!!   ed_cellEdges                   (C) : Holds the edge coordinates for all cells within a block
!!   ed_cellGradNele                (C) : Holds the Nele gradient for all cells within a block
!!   ed_cellGradTele                (C) : Holds the Tele gradient for all cells within a block
!!   ed_cellNele                    (C) : Holds the Nele values for all cells within a block
!!   ed_cellStepTolerance           (R) : Maximum cell fractional error (unit = cell edge) for a ray path step
!!   ed_cellTele                    (C) : Holds the Tele values for all cells within a block
!!   ed_cellVolume                  (C) : Holds the cell's volume for all cells within a block
!!   ed_cellWallThickness           (C) : The cells for each block are treated as having this wall thickness
!!   ed_cellWallThicknessFactor     (R) : This factor x the smallest cell edge gives the final cell thickness
!!   ed_cellZbar                    (C) : Holds the Zbar values for all cells within a block
!!   ed_computeGradNele[X,Y,Z]      (R) : Compute (x,y,z)-components of the number of electron gradients ?
!!   ed_cubicInterpolationZeroDerv  (R) : If true (default), all cubic interpolation vertex derivatives are set to 0
!!   ed_currentStepNumber           (G) : Stores the current step number of the simulation
!!   ed_degrees2rad                 (P) : Angular conversion factor from degrees -> rad
!!   ed_domainErrorMargin[X,Y,Z]    (C) : Computational error margins for domain x,y,z boundaries
!!   ed_domainTolerance             (P) : Relative accuracy for computing the domain error margins
!!   ed_electronCharge              (G) : Electron charge (local copy)
!!   ed_electronMass                (G) : Electron mass (local copy)
!!   ed_energyInTimestep            (C) : Laser energy pumped into the domain at each time step
!!   ed_energyInTotal               (C) : Total laser energy pumped into the domain
!!   ed_energyOutTimestep           (C) : Unused Laser energy exiting the domain at each time step
!!   ed_energyOutTotal              (C) : Total (unused) laser energy exiting the domain
!!   ed_energyProfileFileName       (C) : The file name for writing out the energy profile
!!   ed_energyProfileFileUnit       (C) : The file unit number for the file where the energy profile is written
!!   ed_enforcePositiveNele         (R) : Rescale the Nele gradient such that it is always >= 0?
!!   ed_enforcePositiveTele         (R) : Rescale the Tele gradient such that it is always >= 0?
!!   ed_globalComm                  (C) : Global MPI communicator
!!   ed_globalMe                    (C) : Rank of global processor
!!   ed_globalNumProcs              (C) : Number of processors of the global communicator
!!   ed_gradOrder                   (R) : Number of electrons / electron temperature gradient order
!!   ed_gridGeometry                (C) : Handle to identify the geometrical setup of the calculation
!!   ed_infinitePower               (C) : The infinite power value
!!   ed_infiniteTime                (C) : The infinite time value
!!   ed_infiniteSpeed               (C) : The infinite speed value
!!   ed_Joule2erg                   (P) : Energy conversion factor from Joule -> erg
!!   ed_largestPositiveInteger      (C) : The largest positive integer representable on the current run
!!   ed_largestPositiveReal         (C) : The largest positive real number representable on the current run
!!   ed_laser3Din2D                 (R) : If true, it activates 3D in 2D ray tracing
!!   ed_laser3Din2DwedgeAngle       (R) : The angle (degress, must be < 180) of the wedge for 3D in 2D ray tracing
!!   ed_laser3Din2DwedgeCosine      (C) : The cosine of the 3D in 2D ray tracing rotation angle (same as wedge angle)
!!   ed_laser3Din2DwedgeSine        (C) : The sine of the 3D in 2D ray tracing rotation angle (same as wedge angle)
!!   ed_laser3Din2DwedgeSlope       (C) : The slope m in y = mx, defining the wedge sides for 3D in 2D ray tracing
!!   ed_laserIOMaxNumberOfPositions (R) : Maximum number of positions to store for each IO ray
!!   ed_laserIOMaxNumberOfRays      (R) : Maximum number of IO rays to write out accross each process
!!   ed_laserIONumberOfPositions    (C) : Actual # of positions for each IO ray
!!   ed_laserIONumberOfRaysWritten  (C) : Actual # of IO rays that this process will write out
!!   ed_laserIORayFrequency         (C) : Write out every 'ed_laserIORayFrequency'-th ray
!!   ed_laserIORayPositions         (C) : Stores all x,y,z positions of the IO rays
!!   ed_laserIORayPower             (C) : Stores the power of each of the IO rays
!!   ed_laserIORayTags              (C) : Stores the tags of each of the IO rays
!!   ed_laserIOWrite                (C) : Logical keyword activating IO writeout of particular rays
!!   ed_maxRayCount                 (R) : Maximum number of rays per processor expected
!!   ed_maxPulseSections            (R) : Maximum number of [time|power] sections per pulse
!!   ed_meshComm                    (C) : Mesh MPI communicator
!!   ed_meshMe                      (C) : Rank of mesh processor
!!   ed_meshNumProcs                (C) : Number of processors of the mesh communicator
!!   ed_microns2cm                  (P) : Length conversion factor from microns -> cm
!!   ed_normalizedTolerance         (P) : Maximum allowed deviation from 1.0 for normalization
!!   ed_notSetInteger               (C) : Extreme value (-huge) used for initializing integers
!!   ed_notSetReal                  (C) : Extreme value (-huge) used for initializing reals
!!   ed_numberOfBeams               (R) : Number of beams to set up the laser environment
!!   ed_numberOfPulses              (R) : Number of different kinds of pulses to be used
!!   ed_numberOfSavedRays           (C) : Number of rays saved into the saved rays array
!!   ed_orthogonalTolerance         (P) : Maximum allowed deviation from 0.0 for orthogonalization
!!   ed_particleIndexCount          (C) : Number of particle indices to be used for particle block moving
!!   ed_particleIndexList           (C) : List of particle indices to be used for particle block moving
!!   ed_powerStepTolerance          (R) : Maximum power fractional error (unit = current power) for a ray path step
!!   ed_printBeams                  (R) : Logical keyword activating printing of beam data
!!   ed_printMain                   (R) : Logical keyword activating printing of main data
!!   ed_printPulses                 (R) : Logical keyword activating printing of pulse data
!!   ed_printRays                   (R) : Logical keyword activating printing of ray data
!!   ed_pulses                      (R) : An array of special type, containing all the pulses info (see below)
!!   ed_pulsesAreSetup              (C) : Logical indicator showing if the pulses have been set up
!!   ed_pulseNumberOfSections       (R) : Number of [time|power] sections for each pulse
!!   ed_rayBlockID                  (C) : All different ray block ID's on a processor
!!   ed_rayBlockIDCount             (C) : Number of different ray block ID's on a processor
!!   ed_rayCount                    (C) : Ray count on each processor at any stage of the calculation
!!   ed_rayNumberBlockID            (C) : Number of rays for each ray block ID on a processor
!!   ed_rays                        (C) : Two dimensional array to keep information about each ray (see below)
!!   ed_raysMovedIntoDomain         (C) : Indicates, if the rays are considered to have moved into the domain
!!   ed_raysSaved                   (C) : Array of special type. Will have info about all rays leaving the domain
!!   ed_rayZeroPower                (R) : A ray having power below this value is considered powerless
!!   ed_RungeKuttaMethod            (R) : Specifies the Runge Kutta method to be used for ray tracing
!!   ed_saveOutOfDomainRays         (R) : Should the rays exiting the domain be saved (for diagnostics) ?
!!   ed_speedOfLight                (G) : Speed of light (local copy)
!!   ed_speedOfLightSquared         (C) : The square of the speed of light
!!   ed_tagStart                    (C) : Initial value for tags (always 0, except on 2nd pass split solver)
!!   ed_threadRayTrace              (R) : If true, it activates the threaded section of the ray tracing
!!   ed_unitRoundoff                (C) : The machine epsilon. This + 1.0 is different from 1.0
!!   ed_useEnergyDeposition         (R) : Switch to turn on/off energy deposition computation at run time
!!   ed_useLaserIO                  (R) : Switch to turn on/off writing IO ray data
!!   ed_[x,y,z|min,max]Domain       (R) : The domain bounding box
!!
!!
!!  Explanation of the entries for the rays, pulses and beams
!!  ---------------------------------------------------------
!!
!!   array ed_rays (index1 , index2):
!!
!!     index1                    (C) : Integer handle for the ray properties
!!     index2                    (C) : Ray counting index
!!
!!   type raySaveType:
!!
!!     rayX                      (C) : The global x-coordinate of the saved ray
!!     rayY                      (C) : The global y-coordinate of the saved ray
!!     rayZ                      (C) : The global z-coordinate of the saved ray
!!     rayPower                  (C) : The power of the saved ray
!!     rayTag                    (C) : The global tag of the saved ray
!!
!!   type pulseType:
!!
!!     pulsePower                (R) : The power of the pulse at the pulse time below
!!     pulseTime                 (R) : The pulse time in terms of simulation time
!!
!!   type beamType:
!!
!!     crossSectionFunctionType  (R) : Function type to determine weights in the beams cross section
!!     dimensionality            (C) : Dimensionality of the beam
!!     distanceLens2Target       (C) : Will contain the distance from the lens to the target
!!     frequency                 (C) : Laser frequency of beam in Hz
!!     gaussianCenterMajor       (R) : The gaussian center location along the major semiaxis in cm
!!     gaussianCenterMinor       (R) : The gaussian center location along the minor semiaxis in cm
!!     gaussianExponent          (R) : The gaussian exponent for the beam cross section
!!     gaussianRadiusMajor       (R) : The gaussian radius (e-folding length) along the major semiaxis in cm
!!     gaussianRadiusMinor       (R) : The gaussian radius (e-folding length) along the minor semiaxis in cm
!!     gridDelta1stDim         (R,C) : The tic spacing along the 1st grid dimension
!!     gridDelta2ndDim         (R,C) : The tic spacing along the 2nd grid dimension
!!     gridFirstTic1stDim        (C) : The 1st tic position along the 1st grid dimension
!!     gridFirstTic2ndDim        (C) : The 1st tic position along the 2nd grid dimension
!!     gridnTics1stDim         (R,C) : The # of tic (ray) positions along the 1st grid dimension
!!     gridnTics2ndDim         (R,C) : The # of tic (ray) positions along the 2nd grid dimension
!!     gridSeed                  (C) : Integer seed for random tic (ray) positions
!!     gridSeedMaximum           (C) : The maximum allowed values for the integer seed
!!     gridSeedStepping          (C) : The stepping value for the seed
!!     gridType                  (R) : The type of beam grid to be created for launching the rays
!!     gridWeight                (C) : The beam grid weight (sum of grid point weights)
!!     ignoreBoundaryCondition   (R) : Ignore the domain boundary conditions when creating rays ?
!!     initialRaySpeed           (R) : Initial ray speed in units of light speed
!!     lensSemiAxisMajor         (R) : Length of major lens semiaxis of the beam in cm
!!     lensSemiAxisMinor         (C) : Length of minor lens semiaxis of the beam in cm
!!     lensX                     (R) : Lens center x-coordinate
!!     lensY                     (R) : Lens center y-coordinate
!!     lensZ                     (R) : Lens center z-coordinate
!!     numberOfRays              (R) : Number of rays in beam (can be overwritten in 3D case)
!!     pulseNumber               (R) : Pulse ID number, characterizing beam power in time
!!     pulseStartingTime         (C) : Pulse starting time
!!     pulseEndingTime           (C) : Pulse ending time
!!     semiAxisMajorTorsionAngle (R) : Angle wrt torsion axis of major (target + lens) semiaxis
!!     semiAxisMajorTorsionAxis  (R) : Torsion axis of major (target + lens) semiaxis
!!     semiAxisUnitMajorX        (C) : Unit vector x-component of major (target + lens) semiaxis
!!     semiAxisUnitMajorY        (C) : Unit vector y-component of major (target + lens) semiaxis
!!     semiAxisUnitMajorZ        (C) : Unit vector z-component of major (target + lens) semiaxis
!!     semiAxisUnitMinorX        (C) : Unit vector x-component of minor (target + lens) semiaxis
!!     semiAxisUnitMinorY        (C) : Unit vector y-component of minor (target + lens) semiaxis
!!     semiAxisUnitMinorZ        (C) : Unit vector z-component of minor (target + lens) semiaxis
!!     target2LensMagnification  (C) : The magnification factor for target -> lens data
!!     targetSemiAxisMajor       (R) : Length of major target semiaxis of the beam in cm
!!     targetSemiAxisMinor       (R) : Length of minor target semiaxis of the beam in cm
!!     targetX                   (R) : Target center x-coordinate
!!     targetY                   (R) : Target center y-coordinate
!!     targetZ                   (R) : Target center z-coordinate
!!     wavelength                (R) : Laser wavelength of beam in cm
!!
!!***

Module EnergyDeposition_data

  implicit none

#include "constants.h"
#include "EnergyDeposition.h"
#include "Flash.h"

#ifdef FLASH_GRID_PARTICLES
# include "GridParticles.h"
#endif

  character (len = MAX_STRING_LENGTH), save :: ed_baseName
  character (len = MAX_STRING_LENGTH), save :: ed_energyProfileFileName
  character (len = MAX_STRING_LENGTH), save :: ed_RungeKuttaMethod

  logical, save :: ed_beamsAreSetup       = .false.
  logical, save :: ed_computeGradNeleX
  logical, save :: ed_computeGradNeleY
  logical, save :: ed_computeGradNeleZ
  logical, save :: ed_cubicInterpolationZeroDerv
  logical, save :: ed_enforcePositiveNele
  logical, save :: ed_enforcePositiveTele
  logical, save :: ed_laser3Din2D
  logical, save :: ed_laserIOWrite        = .false.
  logical, save :: ed_printBeams
  logical, save :: ed_printMain
  logical, save :: ed_printPulses
  logical, save :: ed_printRays
  logical, save :: ed_pulsesAreSetup      = .false.
  logical, save :: ed_rayDeterminism
  logical, save :: ed_raysMovedIntoDomain
  logical, save :: ed_saveOutOfDomainRays
  logical, save :: ed_threadRayTrace
  logical, save :: ed_useEnergyDeposition
  logical, save :: ed_useLaserIO

  integer, save :: ed_currentStepNumber
  integer, save :: ed_energyProfileFileUnit
  integer, save :: ed_globalComm
  integer, save :: ed_globalMe
  integer, save :: ed_globalNumProcs
  integer, save :: ed_gradOrder
  integer, save :: ed_gridGeometry
  integer, save :: ed_largestPositiveInteger
  integer, save :: ed_laserIOMaxNumberOfPositions
  integer, save :: ed_laserIOMaxNumberOfRays
  integer, save :: ed_laserIONumberOfRaysWritten
  integer, save :: ed_laserIORayFrequency
  integer, save :: ed_maxRayCount
  integer, save :: ed_meshComm
  integer, save :: ed_meshMe
  integer, save :: ed_meshNumProcs
  integer, save :: ed_notSetInteger
  integer, save :: ed_numberOfBeams
  integer, save :: ed_numberOfPulses
  integer, save :: ed_numberOfSavedRays
  integer, save :: ed_particleIndexCount
  integer, save :: ed_rayBlockIDCount
  integer, save :: ed_rayCount
  integer, save :: ed_tagStart = 0

#ifdef FLASH_GRID_PARTICLES
  integer, save :: ed_particleIndexList (1:GRPT_ALL)
#endif

  real,    save :: ed_Avogadro
  real,    save :: ed_Boltzmann
  real,    save :: ed_cellStepTolerance
  real,    save :: ed_cellWallThickness
  real,    save :: ed_cellWallThicknessFactor
  real,    save :: ed_domainErrorMarginX
  real,    save :: ed_domainErrorMarginY
  real,    save :: ed_domainErrorMarginZ
  real,    save :: ed_electronCharge
  real,    save :: ed_electronMass
  real,    save :: ed_energyInTimestep
  real,    save :: ed_energyInTotal
  real,    save :: ed_energyOutTimestep
  real,    save :: ed_energyOutTotal
  real,    save :: ed_infinitePower
  real,    save :: ed_infiniteTime
  real,    save :: ed_infiniteSpeed
  real,    save :: ed_largestPositiveReal
  real,    save :: ed_laser3Din2DwedgeAngle
  real,    save :: ed_laser3Din2DwedgeCosine
  real,    save :: ed_laser3Din2DwedgeSine
  real,    save :: ed_laser3Din2DwedgeSlope
  real,    save :: ed_notSetReal
  real,    save :: ed_powerStepTolerance
  real,    save :: ed_rayZeroPower
  real,    save :: ed_speedOfLight
  real,    save :: ed_speedOfLightSquared
  real,    save :: ed_unitRoundoff
  real,    save :: ed_xmaxDomain
  real,    save :: ed_xminDomain
  real,    save :: ed_ymaxDomain
  real,    save :: ed_yminDomain
  real,    save :: ed_zmaxDomain
  real,    save :: ed_zminDomain

  real, parameter :: ed_badTorsionAxis      = 1.e-4
  real, parameter :: ed_degrees2rad         = 1.74532925199433e-2
  real, parameter :: ed_domainTolerance     = 1.e-12
  real, parameter :: ed_Joule2erg           = 1.e+7
  real, parameter :: ed_microns2cm          = 1.e-4
  real, parameter :: ed_normalizedTolerance = 1.e-10
  real, parameter :: ed_orthogonalTolerance = 1.e-10

  type beamType
    character (len = BEAM_STRING_LENGTH) :: crossSectionFunctionType
    integer                              :: dimensionality
    real                                 :: distanceLens2Target
    real                                 :: frequency
    real                                 :: gaussianCenterMajor
    real                                 :: gaussianCenterMinor
    real                                 :: gaussianExponent
    real                                 :: gaussianRadiusMajor
    real                                 :: gaussianRadiusMinor
    real                                 :: gridDelta1stDim
    real                                 :: gridDelta2ndDim
    real                                 :: gridFirstTic1stDim
    real                                 :: gridFirstTic2ndDim
    integer                              :: gridnTics1stDim
    integer                              :: gridnTics2ndDim
    integer                              :: gridSeed
    integer                              :: gridSeedMaximum
    integer                              :: gridSeedStepping
    character (len = BEAM_STRING_LENGTH) :: gridType
    real                                 :: gridWeight
    logical                              :: ignoreBoundaryCondition
    real                                 :: initialRaySpeed
    real                                 :: lensSemiAxisMajor
    real                                 :: lensSemiAxisMinor
    real                                 :: lensX
    real                                 :: lensY
    real                                 :: lensZ
    integer                              :: numberOfRays
    integer                              :: pulseNumber
    real                                 :: pulseStartingTime
    real                                 :: pulseEndingTime
    real                                 :: semiAxisMajorTorsionAngle
    character (len = 1)                  :: semiAxisMajorTorsionAxis
    real                                 :: semiAxisUnitMajorX
    real                                 :: semiAxisUnitMajorY
    real                                 :: semiAxisUnitMajorZ
    real                                 :: semiAxisUnitMinorX
    real                                 :: semiAxisUnitMinorY
    real                                 :: semiAxisUnitMinorZ
    real                                 :: target2LensMagnification
    real                                 :: targetSemiAxisMajor
    real                                 :: targetSemiAxisMinor
    real                                 :: targetX
    real                                 :: targetY
    real                                 :: targetZ
    real                                 :: wavelength
  end type beamType

  type pulseType
    real      :: pulsePower
    real      :: pulseTime
  end type pulseType

  type raySaveType
    real      :: rayX
    real      :: rayY
    real      :: rayZ
    real      :: rayPower
    integer   :: rayTag
  end type raySaveType

  type (beamType),     save, allocatable :: ed_beams                    (:)
  integer,             save, allocatable :: ed_laserIONumberOfPositions (:)
  integer,             save, allocatable :: ed_laserIORayTags           (:)
  integer,             save, allocatable :: ed_pulseNumberOfSections    (:)
  integer,             save, allocatable :: ed_rayBlockID               (:)
  integer,             save, allocatable :: ed_rayNumberBlockID         (:)
  type (raySaveType),  save, allocatable :: ed_raysSaved                (:)

  real,                save, allocatable :: ed_cellCenters              (:,:)
  real,                save, allocatable :: ed_cellEdges                (:,:)
  real,                save, allocatable :: ed_laserIORayPower          (:,:)
  type (pulseType),    save, allocatable :: ed_pulses                   (:,:)
  real,                save, allocatable :: ed_rays                     (:,:)

  real,                save, allocatable :: ed_cellDensity              (:,:,:)
  real,                save, allocatable :: ed_cellNele                 (:,:,:)
  real,                save, allocatable :: ed_cellTele                 (:,:,:)
  real,                save, allocatable :: ed_cellVolume               (:,:,:)
  real,                save, allocatable :: ed_cellZbar                 (:,:,:)
  real,                save, allocatable :: ed_laserIORayPositions      (:,:,:)

  real,                save, allocatable :: ed_cellCubicNele            (:,:,:,:)
  real,                save, allocatable :: ed_cellCubicTele            (:,:,:,:)
  real,                save, allocatable :: ed_cellGradNele             (:,:,:,:)
  real,                save, allocatable :: ed_cellGradTele             (:,:,:,:)

  integer,save  :: ed_gcMaskSize=0
  logical,dimension(:),allocatable,save :: ed_gcMask

  logical, save :: ed_depoVarValid = .FALSE.
  logical, save :: ed_depoVarIsPerMass
  integer, save :: ed_savedDepoStepNumber
  integer, save :: ed_depoReuseMaxSteps
  real,    save :: ed_prevDt = 0.0  !should trigger division by zero if inappropriately used
  real,    save :: ed_depoDt = 0.0  !should trigger division by zero if inappropriately used
  integer, save :: ed_depoVar
  integer, save :: ed_irradVar
  character(len=4), save :: ed_depoVarName
  character(len=4), save :: ed_irradVarName

end Module EnergyDeposition_data
