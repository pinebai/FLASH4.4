!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonImaging_data
!!
!! NAME
!!
!!  ProtonImaging_data
!!
!! SYNOPSIS
!!
!!  use ProtonImaging_data
!!  
!! DESCRIPTION
!!
!!  Data module for ProtonImaging
!!  -----------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (C) means data that is (C)alculated internally by the proton imaging code
!!
!!   pi_Avogadro                    (G) : Avogadro's constant (local copy)
!!   pi_badTiltingAxis              (P) : An internal criterion for having chosen a wrong tilting axis
!!   pi_baseName                    (R) : The simulation base name
!!   pi_beams                     (R,C) : An array of special type, containing all the beams info (see below)
!!   pi_beamsAreSetup               (C) : Logical indicator showing if the proton beams have been set up
!!   pi_Boltzmann                   (G) : Boltzmann's constant (local copy)
!!   pi_cellBfield                  (C) : Holds the magnetic field (x,y,z)-components for all cells within a block
!!   pi_cellBoundary                (C) : Holds the boundary indicators for all cells within a block
!!   pi_cellCurlBfield              (C) : Holds the magnetic curl field (x,y,z)-components for all cells within a block
!!   pi_cellEdges[X,Y,Z]            (C) : Holds the edge (x,y,z)-coordinates for all cells within a block
!!   pi_cellEfield                  (C) : Holds the electric field (x,y,z)-components for all cells within a block
!!   pi_cellStepTolerance           (R) : The allowed cell fractional error (units = cell edge) for a proton path step
!!   pi_cellWallThickness           (C) : The cells for each block are treated as having this wall thickness
!!   pi_cellWallThicknessFactor     (R) : This factor x the smallest cell edge gives the final cell wall thickness
!!   pi_degrees2rad                 (P) : Angular conversion factor from degrees -> rad
!!   pi_detectorDGwriteFormat       (R) : Format string for writing out diagnostic variables to detector file(s)
!!   pi_detectorFileNameTimeStamp   (R) : If true (default), a time stamp is added to each detector file name
!!   pi_detectorFilesID             (C) : Holds the detector files ID number
!!   pi_detectorLNwriteFormat       (C) : Format string for writing out a data line to detector file(s)
!!   pi_detectorsAreSetup           (C) : Logical indicator showing if the proton detectors have been set up
!!   pi_detectorXYwriteFormat       (R) : Format string to write proton (x,y) pairs to detector file(s)
!!   pi_domainErrorMargin[X,Y,Z]    (C) : Computational error margins for domain x,y,z boundaries
!!   pi_domainTolerance             (P) : Relative accuracy for computing the domain error margins
!!   pi_flagDomainMissingProtons    (R) : Should domain missing protons be flagged (program aborted)?
!!   pi_globalComm                  (C) : Global MPI communicator
!!   pi_globalMe                    (C) : Rank of global processor
!!   pi_globalNumProcs              (C) : Number of processors of the global communicator
!!   pi_gridGeometry                (C) : Handle to identify the geometrical setup of the calculation
!!   pi_ignoreElectricalField       (R) : Should the electric field be ignored when protons move through domain?
!!   pi_infiniteTime                (C) : The infinite time value
!!   pi_infiniteSpeed               (C) : The infinite speed value
!!   pi_IOaddBeamCapsules           (R) : If true, the frame of the beam capsule(s) will be added to the plot
!!   pi_IOaddDetectorScreens        (R) : If true, the frame of the detector screen(s) will be added to the plot
!!   pi_IOaddProtonsCapsule2Domain  (R) : If true, the proton path from capsule to domain will be added to the plot
!!   pi_IOaddProtonsDomain2Screen   (R) : If true, the proton path from domain to screen will be added to the plot
!!   pi_IOmaxBlockCrossingNumber    (R) : The (estimated) maximum number of complete block crossings for each proton 
!!   pi_IOmaxPointsPerBlock         (C) : The maximum number of points that can be plotted for each proton per block
!!   pi_IOmaxProtonCount            (C) : Maximum number of IO protons expected per processor
!!   pi_IOnumberOfProtons2Plot      (R) : Number of IO protons that are to be plotted
!!   pi_IOplotProtons               (C) : If true, the current simulation step plots the IO protons
!!   pi_IOprotonCount               (C) : IO proton count on each processor at any stage of the calculation
!!   pi_IOprotonPointCount          (C) : The number of points plotted for each IO proton
!!   pi_IOprotonPoints              (C) : The points plotted (currently only x,y,z positions) for each IO proton
!!   pi_IOprotonTags                (C) : The tags of the protons plotted for each IO proton
!!   pi_IOprotonWriteModulo         (C) : Only protons with tags multiple of this modulo base will be plotted
!!   pi_Joule2erg                   (P) : Energy conversion factor from Joule -> erg
!!   pi_largestPositiveInteger      (C) : The largest positive integer representable on the current run
!!   pi_largestPositiveReal         (C) : The largest positive real number representable on the current run
!!   pi_3Din2D                      (R) : If true, it activates 3D in 2D proton tracing
!!   pi_3Din2DwedgeAngle            (R) : The angle (degress, < 180) of the wedge for 3D in 2D proton tracing
!!   pi_3Din2DwedgeCosine           (C) : The cosine of the 3D in 2D proton tracing wedge angle
!!   pi_3Din2DwedgeSine             (C) : The sine of the 3D in 2D proton tracing wedge angle
!!   pi_3Din2DwedgeSlope            (C) : The slope m in y = mx, defining the wedge sides for 3D in 2D proton tracing
!!   pi_maxProtonCount              (R) : Maximum number of protons and screen protons per processor expected
!!   pi_meshComm                    (C) : Mesh MPI communicator
!!   pi_meshMe                      (C) : Rank of mesh processor
!!   pi_meshNumProcs                (C) : Number of processors of the mesh communicator
!!   pi_MeV2erg                     (P) : Energy conversion factor from MeV -> erg
!!   pi_microns2cm                  (P) : Length conversion factor from microns -> cm
!!   pi_mpiScreenProtonType         (C) : Screen proton type structure for mpi communications
!!   pi_normalizedTolerance         (P) : Maximum allowed deviation from 1.0 for normalization on unit vectors
!!   pi_notSetInteger               (C) : Extreme value (-huge) used for initializing integers
!!   pi_notSetReal                  (C) : Extreme value (-huge) used for initializing reals
!!   pi_numberOfBeams               (R) : Number of beams to set up the proton imaging environment
!!   pi_numberOfDetectors           (R) : Number of detectors to set up the proton imaging environment
!!   pi_numberOfDiagnostics         (P) : Number of diagnostic variables currently evaluated by the code
!!   pi_opaqueBoundaries            (R) : If true, the protons do not go through cells belonging to boundaries
!!   pi_orthogonalTolerance         (P) : Maximum allowed deviation from 0.0 for orthogonalization on unit vectors
!!   pi_particleIndexCount          (C) : Number of particle indices to be used for particle block moving
!!   pi_particleIndexList           (C) : List of particle indices to be used for particle block moving
!!   pi_printBeams                  (R) : Logical keyword activating printing of beam data
!!   pi_printDetectors              (R) : Logical keyword activating printing of detector data
!!   pi_printMain                   (R) : Logical keyword activating printing of main data
!!   pi_printProtons                (R) : Logical keyword activating printing of proton data
!!   pi_protonBlockID               (C) : All different proton block ID's on a processor
!!   pi_protonBlockIDCount          (C) : Number of different proton block ID's on a processor
!!   pi_protonCharge                (G) : Proton charge (local copy)
!!   pi_protonChargePerMass         (C) : Proton charge per mass
!!   pi_protonCount                 (C) : Proton count on each processor at any stage of the calculation
!!   pi_protonDeterminism           (R) : If true, Grid Unit will use the Sieve Algorithm to move proton particles
!!   pi_protonMass                  (G) : Proton mass (local copy)
!!   pi_protonNumberBlockID         (C) : Number of protons for each proton block ID on a processor
!!   pi_protons                     (C) : Two dimensional array to keep information about each proton (see below)
!!   pi_protonsMovedIntoDomain      (C) : Indicates, if the protons are considered to have moved into the domain
!!   pi_randomNumberSeedArray       (C) : Will hold the seeds for the random number generator
!!   pi_randomNumberSeedIncrement   (R) : Sets the seed increment for the random number generator
!!   pi_randomNumberSeedInitial     (R) : Sets the initial seeds for the random number generator
!!   pi_recalculateCellData         (R) : If true, the proton imaging calculates its own cell data for each block
!!   pi_recordOffScreenProtons      (R) : If true, the protons missing the detector screen will also be recorded.
!!   pi_RungeKuttaMethod            (R) : Specifies the Runge Kutta method to be used for proton tracing.
!!   pi_screenProtonBucketCount     (C) : Counts how many screen protons are currently in screen proton bucket
!!   pi_screenProtonBucketSize      (R) : Bucket size for flushing out screen protons to disk
!!   pi_screenProtonCount           (C) : Screen proton count on each processor at any stage of the calculation
!!   pi_screenProtonCountOffsets    (C) : Screen proton count offsets info for master processor
!!   pi_screenProtonCountProcs      (C) : Screen proton count for each processor info for master processor
!!   pi_screenProtonDiagnostics     (R) : If true, calculates/records extra diagnostic values for the screen protons
!!   pi_screenProtonRecordCount     (C) : Screen proton record count (currently 6: x,y,Jx,Kx,Ky,Kz) 
!!   pi_screenProtons               (C) : Screen proton array to record protons on detector screen
!!   pi_speedOfLight                (G) : Speed of light (local copy)
!!   pi_speedOfLightInv             (C) : The inverse of the speed of light
!!   pi_speedOfLightSquared         (C) : The square of the speed of light
!!   pi_squareRoot4Pi               (C) : The square root of 4 pi (to convert FLASH magnetic B to Gauss units)
!!   pi_tagMax                      (C) : The current maximum proton tag value at any stage of the calculation
!!   pi_threadProtonTrace           (R) : If true, it activates the threaded section of the proton tracing
!!   pi_totalProtons2Launch         (C) : The total number of protons to be launched (over all beams)
!!   pi_unitRoundoff                (C) : The machine epsilon. This + 1.0 is different from 1.0
!!   pi_useIOprotonPlot             (R) : If true, proton data will be written to plot files
!!   pi_useParabolicApproximation   (R) : If true, the parabolic path approximation is used whenever possible
!!   pi_useProtonImaging            (R) : Switch to turn on/off proton imaging at run time
!!   pi_[x,y]Circle                 (C) : Will contain (x,y) pairs in a circle (beam target areas)
!!   pi_[x,y,z|min,max]Domain       (R) : The domain bounding box
!!   pi_[x,y,z]Sphere               (C) : Will contain (x,y,z) triples in a sphere (beam capsules)
!!
!!
!!  Explanation of the entries for the protons, screen protons, beams and detectors
!!  -------------------------------------------------------------------------------
!!
!!   array pi_protons (index1 , index2):
!!
!!     index1                    (C) : Integer handle for the proton properties
!!     index2                    (C) : Proton counting index
!!
!!   array pi_screenProtons (index1 , index2):
!!
!!     index1                    (C) : Integer handle for the screen proton properties
!!     index2                    (C) : Screen proton counting index
!!
!!   type beamType:
!!
!!     apertureAngle             (R) : Aperture angle of the conical beam
!!     axisUnit1X                (C) : 1st axis unit vector x-component
!!     axisUnit1Y                (C) : 1st axis unit vector y-component
!!     axisUnit1Z                (C) : 1st axis unit vector z-component
!!     axisUnit2X                (C) : 2nd axis unit vector x-component
!!     axisUnit2Y                (C) : 2nd axis unit vector y-component
!!     axisUnit2Z                (C) : 2nd axis unit vector z-component
!!     axisUnit3X                (C) : 3rd axis unit vector x-component
!!     axisUnit3Y                (C) : 3rd axis unit vector y-component
!!     axisUnit3Z                (C) : 3rd axis unit vector z-component
!!     capsuleGrainIndexI        (C) : The capsule grain i-index of the index triple (i,j,k)
!!     capsuleGrainIndexJ        (C) : The capsule grain j-index of the index triple (i,j,k)
!!     capsuleGrainIndexK        (C) : The capsule grain k-index of the index triple (i,j,k)
!!     capsuleGrainLevel         (R) : The capsule grain level
!!     capsuleGrainLocalX        (C) : The capsule grain local x-coordinate
!!     capsuleGrainLocalY        (C) : The capsule grain local y-coordinate
!!     capsuleGrainLocalZ        (C) : The capsule grain local z-coordinate
!!     capsuleGrainProtons       (C) : The capsule grain number of protons
!!     capsuleGrainSize          (C) : The size of each capsule grain (currently the side length of a cube)
!!     capsuleNumberOfGrains     (C) : The capsule total number of grains
!!     capsuleRadius             (R) : Capsule radius in cm
!!     capsuleX                  (R) : Capsule center x-coordinate
!!     capsuleY                  (R) : Capsule center y-coordinate
!!     capsuleZ                  (R) : Capsule center z-coordinate
!!     detector                  (R) : Target detector nr of the beam
!!     dimensionality            (C) : Dimensionality of the beam
!!     distanceCapsule2Target    (C) : Will contain the distance from the capsule center to the target
!!     initialProtonSpeed        (C) : Initial proton speed in cm/s
!!     noBoundaryCondition       (R) : Ignore the domain boundary conditions when creating protons ?
!!     numberOfProtons         (R,C) : Desired number of protons in beam (may be adjusted by the code)
!!     numberOfProtonsPerGrain   (C) : Number of protons per capsule grain in beam (fixed througout)
!!     protonEnergy              (R) : Energy of protons in beam in MeV
!!     randomNumberSeed          (C) : Integer seed for random number generator
!!     targetRadius            (R,C) : Length of target radius (calculated from aperture angle, if not set)
!!     targetX                   (R) : Target center x-coordinate
!!     targetY                   (R) : Target center y-coordinate
!!     targetZ                   (R) : Target center z-coordinate
!!     time2Launch               (R) : If the simulation time >= time2Launch -> the beam launches its protons
!!
!!   type detectorType:
!!
!!     alignWRTbeamNr            (R) : place detector screen along what beam nr (if <= 0, no placing)
!!     axisXunitX                (C) : detector X axis unit vector x-component
!!     axisXunitY                (C) : detector X axis unit vector y-component
!!     axisXunitZ                (C) : detector X axis unit vector z-component
!!     axisYunitX                (C) : detector Y axis unit vector x-component
!!     axisYunitY                (C) : detector Y axis unit vector y-component
!!     axisYunitZ                (C) : detector Y axis unit vector z-component
!!     centerX                 (R,C) : detector center x-coordinate (calculated, if aligned wrt beam)
!!     centerY                 (R,C) : detector center y-coordinate (calculated, if aligned wrt beam)
!!     centerZ                 (R,C) : detector center z-coordinate (calculated, if aligned wrt beam)
!!     cornerLowerLeftX          (C) : x-coordinate lower  left corner ([-1,-1] direction axis X,Y unit vectors)
!!     cornerLowerLeftY          (C) : y-coordinate lower  left corner ([-1,-1] direction axis X,Y unit vectors)
!!     cornerLowerLeftZ          (C) : z-coordinate lower  left corner ([-1,-1] direction axis X,Y unit vectors)
!!     cornerLowerRightX         (C) : x-coordinate lower right corner ([+1,-1] direction axis X,Y unit vectors)
!!     cornerLowerRightY         (C) : y-coordinate lower right corner ([+1,-1] direction axis X,Y unit vectors)
!!     cornerLowerRightZ         (C) : z-coordinate lower right corner ([+1,-1] direction axis X,Y unit vectors)
!!     cornerUpperLeftX          (C) : x-coordinate upper  left corner ([-1,+1] direction axis X,Y unit vectors)
!!     cornerUpperLeftY          (C) : y-coordinate upper  left corner ([-1,+1] direction axis X,Y unit vectors)
!!     cornerUpperLeftZ          (C) : z-coordinate upper  left corner ([-1,+1] direction axis X,Y unit vectors)
!!     cornerUpperRightX         (C) : x-coordinate upper right corner ([+1,+1] direction axis X,Y unit vectors)
!!     cornerUpperRightY         (C) : y-coordinate upper right corner ([+1,+1] direction axis X,Y unit vectors)
!!     cornerUpperRightZ         (C) : z-coordinate upper right corner ([+1,+1] direction axis X,Y unit vectors)
!!     distance2beamCapsule      (R) : distance from detector center to beam capsule center (if aligned wrt beam)
!!     normalX                 (R,C) : detector normal vector x-coordinate (calculated, if aligned wrt beam)
!!     normalY                 (R,C) : detector normal vector y-coordinate (calculated, if aligned wrt beam)
!!     normalZ                 (R,C) : detector normal vector z-coordinate (calculated, if aligned wrt beam)
!!     pinholeDist2Det           (R) : pinhole center distance from the detector center
!!     pinholeRadius             (R) : pinhole radius
!!     pinholeX                  (C) : pinhole center x-coordinate
!!     pinholeY                  (C) : pinhole center y-coordinate
!!     pinholeZ                  (C) : pinhole center z-coordinate
!!     sideLength                (R) : side length of the square detector (in cm)
!!     sideLengthHalf            (C) : half of the side length of the square detector
!!     sideLengthInv             (C) : inverse of the side length of the square detector
!!     sideTiltingAngle          (R) : side tilting angle from tilting axis (in degrees)
!!     sideTiltingAxis           (R) : global coordinate axis to be used for side tilting ('x','y' or 'z')
!!
!!***

Module ProtonImaging_data

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

#ifdef FLASH_GRID_PARTICLES
# include "GridParticles.h"
#endif

  character (len = MAX_STRING_LENGTH), save :: pi_baseName
  character (len = 10               ), save :: pi_detectorDGwriteFormat
  character (len = 24               ), save :: pi_detectorLNwriteFormat
  character (len = 10               ), save :: pi_detectorXYwriteFormat
  character (len = MAX_STRING_LENGTH), save :: pi_RungeKuttaMethod

  logical, save :: pi_3Din2D
  logical, save :: pi_beamsAreSetup     = .false.
  logical, save :: pi_detectorFileNameTimeStamp
  logical, save :: pi_detectorsAreSetup = .false.
  logical, save :: pi_flagDomainMissingProtons
  logical, save :: pi_ignoreElectricalField
  logical, save :: pi_IOaddBeamCapsules
  logical, save :: pi_IOaddDetectorScreens
  logical, save :: pi_IOaddProtonsCapsule2Domain
  logical, save :: pi_IOaddProtonsDomain2Screen
  logical, save :: pi_IOplotProtons
  logical, save :: pi_opaqueBoundaries
  logical, save :: pi_printBeams
  logical, save :: pi_printDetectors
  logical, save :: pi_printMain
  logical, save :: pi_printProtons
  logical, save :: pi_protonDeterminism
  logical, save :: pi_protonsMovedIntoDomain
  logical, save :: pi_recalculateCellData
  logical, save :: pi_recordOffScreenProtons
  logical, save :: pi_screenProtonDiagnostics
  logical, save :: pi_threadProtonTrace
  logical, save :: pi_useIOprotonPlot
  logical, save :: pi_useParabolicApproximation
  logical, save :: pi_useProtonImaging

  integer, save :: pi_currentStepNumber
  integer, save :: pi_globalComm
  integer, save :: pi_globalMe
  integer, save :: pi_globalNumProcs
  integer, save :: pi_gridGeometry
  integer, save :: pi_IOmaxBlockCrossingNumber
  integer, save :: pi_IOmaxPointsPerBlock
  integer, save :: pi_IOmaxProtonCount
  integer, save :: pi_IOnumberOfProtons2Plot
  integer, save :: pi_IOprotonCount
  integer, save :: pi_IOprotonWriteModulo
  integer, save :: pi_largestPositiveInteger
  integer, save :: pi_maxProtonCount
  integer, save :: pi_meshComm
  integer, save :: pi_meshMe
  integer, save :: pi_meshNumProcs
  integer, save :: pi_mpiScreenProtonType
  integer, save :: pi_notSetInteger
  integer, save :: pi_numberOfBeams
  integer, save :: pi_numberOfDetectors
  integer, save :: pi_randomNumberSeedIncrement
  integer, save :: pi_randomNumberSeedInitial
  integer, save :: pi_particleIndexCount
  integer, save :: pi_protonBlockIDCount
  integer, save :: pi_protonCount
  integer, save :: pi_screenProtonBucketSize
  integer, save :: pi_screenProtonCount
  integer, save :: pi_screenProtonRecordCount
  integer, save :: pi_tagMax
  integer, save :: pi_totalProtons2Launch

#ifdef FLASH_GRID_PARTICLES
  integer, save :: pi_particleIndexList (1:GRPT_ALL)
#endif

  real,    save :: pi_Avogadro
  real,    save :: pi_Boltzmann
  real,    save :: pi_cellStepTolerance
  real,    save :: pi_cellWallThickness
  real,    save :: pi_cellWallThicknessFactor
  real,    save :: pi_domainErrorMarginX
  real,    save :: pi_domainErrorMarginY
  real,    save :: pi_domainErrorMarginZ
  real,    save :: pi_protonCharge
  real,    save :: pi_protonChargePerMass
  real,    save :: pi_protonMass
  real,    save :: pi_infiniteTime
  real,    save :: pi_infiniteSpeed
  real,    save :: pi_largestPositiveReal
  real,    save :: pi_3Din2DwedgeAngle
  real,    save :: pi_3Din2DwedgeCosine
  real,    save :: pi_3Din2DwedgeSine
  real,    save :: pi_3Din2DwedgeSlope
  real,    save :: pi_notSetReal
  real,    save :: pi_speedOfLight
  real,    save :: pi_speedOfLightInv
  real,    save :: pi_speedOfLightSquared
  real,    save :: pi_squareRoot4Pi
  real,    save :: pi_unitRoundoff
  real,    save :: pi_xmaxDomain
  real,    save :: pi_xminDomain
  real,    save :: pi_ymaxDomain
  real,    save :: pi_yminDomain
  real,    save :: pi_zmaxDomain
  real,    save :: pi_zminDomain

  integer, parameter :: pi_numberOfDiagnostics = 4

  real,    parameter :: pi_badTiltingAxis      = 1.e-4
  real,    parameter :: pi_degrees2rad         = 1.74532925199433e-2
  real,    parameter :: pi_domainTolerance     = 1.e-12
  real,    parameter :: pi_Joule2erg           = 1.e+7
  real,    parameter :: pi_MeV2erg             = 1.602176565e-6
  real,    parameter :: pi_microns2cm          = 1.e-4
  real,    parameter :: pi_normalizedTolerance = 1.e-10
  real,    parameter :: pi_orthogonalTolerance = 1.e-10

  type beamType
    real                :: apertureAngle
    real                :: axisUnit1X
    real                :: axisUnit1Y
    real                :: axisUnit1Z
    real                :: axisUnit2X
    real                :: axisUnit2Y
    real                :: axisUnit2Z
    real                :: axisUnit3X
    real                :: axisUnit3Y
    real                :: axisUnit3Z
    integer             :: capsuleGrainIndexI
    integer             :: capsuleGrainIndexJ
    integer             :: capsuleGrainIndexK
    integer             :: capsuleGrainLevel
    real                :: capsuleGrainLocalX
    real                :: capsuleGrainLocalY
    real                :: capsuleGrainLocalZ
    integer             :: capsuleGrainProtons
    real                :: capsuleGrainSize
    integer             :: capsuleNumberOfGrains
    real                :: capsuleRadius
    real                :: capsuleX
    real                :: capsuleY
    real                :: capsuleZ
    integer             :: detector
    integer             :: dimensionality
    real                :: distanceCapsule2Target
    real                :: initialProtonSpeed
    logical             :: noBoundaryCondition
    integer             :: numberOfProtons
    integer             :: numberOfProtonsPerGrain
    real                :: protonEnergy
    integer             :: randomNumberSeed
    real                :: targetRadius
    real                :: targetX
    real                :: targetY
    real                :: targetZ
    real                :: time2Launch
  end type beamType

  type detectorType
    integer             :: alignWRTbeamNr
    real                :: axisXunitX
    real                :: axisXunitY
    real                :: axisXunitZ
    real                :: axisYunitX
    real                :: axisYunitY
    real                :: axisYunitZ
    real                :: centerX
    real                :: centerY
    real                :: centerZ
    real                :: cornerLowerLeftX
    real                :: cornerLowerLeftY
    real                :: cornerLowerLeftZ
    real                :: cornerLowerRightX
    real                :: cornerLowerRightY
    real                :: cornerLowerRightZ
    real                :: cornerUpperLeftX
    real                :: cornerUpperLeftY
    real                :: cornerUpperLeftZ
    real                :: cornerUpperRightX
    real                :: cornerUpperRightY
    real                :: cornerUpperRightZ
    real                :: distance2beamCapsule
    real                :: normalX
    real                :: normalY
    real                :: normalZ
    real                :: pinholeDist2Det
    real                :: pinholeRadius
    real                :: pinholeX
    real                :: pinholeY
    real                :: pinholeZ
    real                :: sideLength
    real                :: sideLengthHalf
    real                :: sideLengthInv
    real                :: sideTiltingAngle
    character (len = 1) :: sideTiltingAxis
  end type detectorType

  type (beamType),         save, allocatable :: pi_beams                    (:)
  real,                    save, allocatable :: pi_cellEdgesX               (:)
  real,                    save, allocatable :: pi_cellEdgesY               (:)
  real,                    save, allocatable :: pi_cellEdgesZ               (:)
  integer,                 save, allocatable :: pi_detectorFilesID          (:)
  type (detectorType),     save, allocatable :: pi_detectors                (:)
  integer,                 save, allocatable :: pi_IOprotonPointCount       (:)
  integer,                 save, allocatable :: pi_IOprotonTags             (:)
  integer,                 save, allocatable :: pi_protonBlockID            (:)
  integer,                 save, allocatable :: pi_protonNumberBlockID      (:)
  integer,                 save, allocatable :: pi_randomNumberSeedArray    (:)
  integer,                 save, allocatable :: pi_screenProtonBucketCount  (:)   ! only for master processor
  integer,                 save, allocatable :: pi_screenProtonCountOffsets (:)   ! only for master processor
  integer,                 save, allocatable :: pi_screenProtonCountProcs   (:)   ! only for master processor
  real,                    save, allocatable :: pi_xCircle                  (:)
  real,                    save, allocatable :: pi_yCircle                  (:)
  real,                    save, allocatable :: pi_xSphere                  (:)
  real,                    save, allocatable :: pi_ySphere                  (:)
  real,                    save, allocatable :: pi_zSphere                  (:)
  real,                    save, allocatable :: pi_protons                  (:,:)
  real,                    save, allocatable :: pi_screenProtons            (:,:)
  real,                    save, allocatable :: pi_cellBoundary             (:,:,:)
  real,                    save, allocatable :: pi_IOprotonPoints           (:,:,:)
  real,                    save, allocatable :: pi_screenProtonBuckets      (:,:,:) ! (properties, count, detector)
  real,                    save, allocatable :: pi_cellBfield               (:,:,:,:)
  real,                    save, allocatable :: pi_cellCurlBfield           (:,:,:,:)
  real,                    save, allocatable :: pi_cellEfield               (:,:,:,:)

end Module ProtonImaging_data
