!!****ih* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_interface
!!
!! NAME
!!
!!  ed_interface
!!
!! SYNOPSIS
!!
!!   use ed_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Energy Deposition unit that defines its
!!  private interfaces.
!!
!!***

Module ed_interface

  interface
     subroutine ed_beam2DGridPointsRegular (semiAxis,                      &
                                            nTics,                         &
                                            delta,                         &
                                            firstTic,                      &
                                            startGrid,                     &
                                            maxGridPoints,                 &
                                                           moreGridPoints, &
                                                           nGridPoints,    &
                                                           xGrid           )

       real,    intent (in)    :: semiAxis
       integer, intent (in)    :: nTics
       real,    intent (in)    :: delta
       real,    intent (in)    :: firstTic
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
     end subroutine ed_beam2DGridPointsRegular
  end interface

  interface
     subroutine ed_beam2DGridPointsStatistical (semiAxis,                      &
                                                seed,                          &
                                                targetTotalGridPoints,         &
                                                startGrid,                     &
                                                maxGridPoints,                 &
                                                               moreGridPoints, &
                                                               nGridPoints,    &
                                                               xGrid           )

       real,    intent (in)    :: semiAxis
       integer, intent (in)    :: seed
       integer, intent (in)    :: targetTotalGridPoints
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
     end subroutine ed_beam2DGridPointsStatistical
  end interface

  interface
     subroutine ed_beam2DGridSetupRegular (semiAxis,             &
                                           nGridPoints,          &
                                                        nTics,   &
                                                        delta,   &
                                                        firstTic )

       real,    intent (in)   :: semiAxis
       integer, intent (in)   :: nGridPoints
       integer, intent (out)  :: nTics
       real,    intent (out)  :: delta
       real,    intent (out)  :: firstTic
     end subroutine ed_beam2DGridSetupRegular
  end interface

  interface
     subroutine ed_beam2DGridSetupStatistical (seedInitialize,               &
                                               seedIncrement,                &
                                                               seedMaximum,  &
                                                               seedStepping, &
                                                               seed          )

       logical, intent (in)    :: seedInitialize
       logical, intent (in)    :: seedIncrement
       integer, intent (inout) :: seedMaximum
       integer, intent (inout) :: seedStepping
       integer, intent (out)   :: seed
     end subroutine ed_beam2DGridSetupStatistical
  end interface

  interface
     subroutine ed_beam2DGridWeightRegular (gridRadialOrigin,                  &
                                            ux,                                &
                                            target2LensMagnify,                &
                                            semiAxis,                          &
                                            crossSectionFunctionType,          &
                                            gaussianExponent,                  &
                                            gaussianRadius,                    &
                                            gaussianCenter,                    &
                                            nTics,                             &
                                            delta,                             &
                                            firstTic,                          &
                                            totalGridPoints,                   &
                                                                    gridWeight )

       real,              intent (in)  :: gridRadialOrigin
       real,              intent (in)  :: ux
       real,              intent (in)  :: target2LensMagnify
       real,              intent (in)  :: semiAxis
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadius
       real,              intent (in)  :: gaussianCenter
       integer,           intent (in)  :: nTics
       real,              intent (in)  :: delta
       real,              intent (in)  :: firstTic
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam2DGridWeightRegular
  end interface

  interface
     subroutine ed_beam2DGridWeightStatistical (gridRadialOrigin,                  &
                                                semiAxis,                          &
                                                crossSectionFunctionType,          &
                                                gaussianExponent,                  &
                                                gaussianRadius,                    &
                                                gaussianCenter,                    &
                                                seed,                              &
                                                totalGridPoints,                   &
                                                                        gridWeight )

       real,              intent (in)  :: gridRadialOrigin
       real,              intent (in)  :: semiAxis
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadius
       real,              intent (in)  :: gaussianCenter
       integer,           intent (in)  :: seed
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam2DGridWeightStatistical
  end interface

  interface
     subroutine ed_beam3DGridPointsDelta (semiAxisMajor,                      &
                                          semiAxisMinor,                      &
                                          nTicsSemiAxisMajor,                 &
                                          nTicsSemiAxisMinor,                 &
                                          deltaSemiAxisMajor,                 &
                                          deltaSemiAxisMinor,                 &
                                          firstTicSemiAxisMajor,              &
                                          firstTicSemiAxisMinor,              &
                                          startGrid,                          &
                                          maxGridPoints,                      &
                                                              moreGridPoints, &
                                                              nGridPoints,    &
                                                              xGrid,          &
                                                              yGrid           )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (in)    :: nTicsSemiAxisMajor
       integer, intent (in)    :: nTicsSemiAxisMinor
       real,    intent (in)    :: deltaSemiAxisMajor
       real,    intent (in)    :: deltaSemiAxisMinor
       real,    intent (in)    :: firstTicSemiAxisMajor
       real,    intent (in)    :: firstTicSemiAxisMinor
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
       real,    intent (out)   :: yGrid (1:maxGridPoints)
     end subroutine ed_beam3DGridPointsDelta
  end interface

  interface
     subroutine ed_beam3DGridPointsRadial (semiAxisMajor,                   &
                                           semiAxisMinor,                   &
                                           nTicsRadial,                     &
                                           nTicsAngular,                    &
                                           deltaRadial,                     &
                                           deltaAngular,                    &
                                           firstTicRadial,                  &
                                           firstTicAngular,                 &
                                           startGrid,                       &
                                           maxGridPoints,                   &
                                                            moreGridPoints, &
                                                            nGridPoints,    &
                                                            xGrid,          &
                                                            yGrid           )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (in)    :: nTicsRadial
       integer, intent (in)    :: nTicsAngular
       real,    intent (in)    :: deltaRadial
       real,    intent (in)    :: deltaAngular
       real,    intent (in)    :: firstTicRadial
       real,    intent (in)    :: firstTicAngular
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
       real,    intent (out)   :: yGrid (1:maxGridPoints)
     end subroutine ed_beam3DGridPointsRadial
  end interface

  interface
     subroutine ed_beam3DGridPointsRecBeam (nTicsSemiAxisMajor,                 &
                                            nTicsSemiAxisMinor,                 &
                                            deltaSemiAxisMajor,                 &
                                            deltaSemiAxisMinor,                 &
                                            startGrid,                          &
                                            maxGridPoints,                      &
                                                                moreGridPoints, &
                                                                nGridPoints,    &
                                                                xGrid,          &
                                                                yGrid           )

       integer, intent (in)    :: nTicsSemiAxisMajor
       integer, intent (in)    :: nTicsSemiAxisMinor
       real,    intent (in)    :: deltaSemiAxisMajor
       real,    intent (in)    :: deltaSemiAxisMinor
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
       real,    intent (out)   :: yGrid (1:maxGridPoints)
     end subroutine ed_beam3DGridPointsRecBeam
  end interface

  interface
     subroutine ed_beam3DGridPointsSquare (semiAxisMajor,                      &
                                           semiAxisMinor,                      &
                                           nTicsSemiAxisMajor,                 &
                                           nTicsSemiAxisMinor,                 &
                                           delta,                              &
                                           startGrid,                          &
                                           maxGridPoints,                      &
                                                               moreGridPoints, &
                                                               nGridPoints,    &
                                                               xGrid,          &
                                                               yGrid           )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (in)    :: nTicsSemiAxisMajor
       integer, intent (in)    :: nTicsSemiAxisMinor
       real,    intent (in)    :: delta
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
       real,    intent (out)   :: yGrid (1:maxGridPoints)
     end subroutine ed_beam3DGridPointsSquare
  end interface

  interface
     subroutine ed_beam3DGridPointsStatistical (semiAxisMajor,                      &
                                                semiAxisMinor,                      &
                                                seed,                               &
                                                targetTotalGridPoints,              &
                                                startGrid,                          &
                                                maxGridPoints,                      &
                                                                    moreGridPoints, &
                                                                    nGridPoints,    &
                                                                    xGrid,          &
                                                                    yGrid           )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (in)    :: seed
       integer, intent (in)    :: targetTotalGridPoints
       logical, intent (inout) :: startGrid
       integer, intent (in)    :: maxGridPoints
       logical, intent (out)   :: moreGridPoints
       integer, intent (out)   :: nGridPoints
       real,    intent (out)   :: xGrid (1:maxGridPoints)
       real,    intent (out)   :: yGrid (1:maxGridPoints)
     end subroutine ed_beam3DGridPointsStatistical
  end interface

  interface
     subroutine ed_beam3DGridSetupDelta (semiAxisMajor,                             &
                                         semiAxisMinor,                             &
                                         deltaSemiAxisMajor,                        &
                                         deltaSemiAxisMinor,                        &
                                                             nGridPoints,           &
                                                             nTicsSemiAxisMajor,    &
                                                             nTicsSemiAxisMinor,    &
                                                             firstTicSemiAxisMajor, &
                                                             firstTicSemiAxisMinor  )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       real,    intent (in)    :: deltaSemiAxisMajor
       real,    intent (in)    :: deltaSemiAxisMinor
       integer, intent (inout) :: nGridPoints
       integer, intent (out)   :: nTicsSemiAxisMajor
       integer, intent (out)   :: nTicsSemiAxisMinor
       real,    intent (out)   :: firstTicSemiAxisMajor
       real,    intent (out)   :: firstTicSemiAxisMinor
     end subroutine ed_beam3DGridSetupDelta
  end interface

  interface
     subroutine ed_beam3DGridSetupRadial (nGridPoints,    &
                                          nTicsRadial,    &
                                          nTicsAngular,   &
                                          deltaRadial,    &
                                          deltaAngular,   &
                                          firstTicRadial, &
                                          firstTicAngular )

       integer, intent (inout) :: nGridPoints
       integer, intent (inout) :: nTicsRadial
       integer, intent (inout) :: nTicsAngular
       real,    intent (out)   :: deltaRadial
       real,    intent (out)   :: deltaAngular
       real,    intent (out)   :: firstTicRadial
       real,    intent (out)   :: firstTicAngular
     end subroutine ed_beam3DGridSetupRadial
  end interface

  interface
     subroutine ed_beam3DGridSetupRecBeam (semiAxisMajor,      &
                                           semiAxisMinor,      &
                                           nGridPoints,        &
                                           nTicsSemiAxisMajor, &
                                           nTicsSemiAxisMinor, &
                                           deltaSemiAxisMajor, &
                                           deltaSemiAxisMinor  )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (inout) :: nGridPoints
       integer, intent (inout) :: nTicsSemiAxisMajor
       integer, intent (inout) :: nTicsSemiAxisMinor
       real,    intent (out)   :: deltaSemiAxisMajor
       real,    intent (out)   :: deltaSemiAxisMinor
     end subroutine ed_beam3DGridSetupRecBeam
  end interface

  interface
     subroutine ed_beam3DGridSetupSquare (semiAxisMajor,                     &
                                          semiAxisMinor,                     &
                                                         nGridPoints,        &
                                                         nTicsSemiAxisMajor, &
                                                         nTicsSemiAxisMinor, &
                                                         delta               )

       real,    intent (in)    :: semiAxisMajor
       real,    intent (in)    :: semiAxisMinor
       integer, intent (inout) :: nGridPoints
       integer, intent (out)   :: nTicsSemiAxisMajor
       integer, intent (out)   :: nTicsSemiAxisMinor
       real,    intent (out)   :: delta
     end subroutine ed_beam3DGridSetupSquare
  end interface

  interface
     subroutine ed_beam3DGridSetupStatistical (seedInitialize,               &
                                               seedIncrement,                &
                                                               seedMaximum,  &
                                                               seedStepping, &
                                                               seed          )

       logical, intent (in)    :: seedInitialize
       logical, intent (in)    :: seedIncrement
       integer, intent (inout) :: seedMaximum
       integer, intent (inout) :: seedStepping
       integer, intent (out)   :: seed
     end subroutine ed_beam3DGridSetupStatistical
  end interface

  interface
     subroutine ed_beam3DGridWeightDelta (semiAxisMajor,                    &
                                          semiAxisMinor,                    &
                                          crossSectionFunctionType,         &
                                          gaussianExponent,                 &
                                          gaussianRadiusMajor,              &
                                          gaussianRadiusMinor,              &
                                          gaussianCenterMajor,              &
                                          gaussianCenterMinor,              &
                                          nTicsSemiAxisMajor,               &
                                          nTicsSemiAxisMinor,               &
                                          deltaSemiAxisMajor,               &
                                          deltaSemiAxisMinor,               &
                                          firstTicSemiAxisMajor,            &
                                          firstTicSemiAxisMinor,            &
                                          totalGridPoints,                  &
                                                                 gridWeight )

       real,              intent (in)  :: semiAxisMajor
       real,              intent (in)  :: semiAxisMinor
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadiusMajor
       real,              intent (in)  :: gaussianRadiusMinor
       real,              intent (in)  :: gaussianCenterMajor
       real,              intent (in)  :: gaussianCenterMinor
       integer,           intent (in)  :: nTicsSemiAxisMajor
       integer,           intent (in)  :: nTicsSemiAxisMinor
       real,              intent (in)  :: deltaSemiAxisMajor
       real,              intent (in)  :: deltaSemiAxisMinor
       real,              intent (in)  :: firstTicSemiAxisMajor
       real,              intent (in)  :: firstTicSemiAxisMinor
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam3DGridWeightDelta
  end interface

  interface
     subroutine ed_beam3DGridWeightRadial (semiAxisMajor,                    &
                                           semiAxisMinor,                    &
                                           crossSectionFunctionType,         &
                                           gaussianExponent,                 &
                                           gaussianRadiusMajor,              &
                                           gaussianRadiusMinor,              &
                                           gaussianCenterMajor,              &
                                           gaussianCenterMinor,              &
                                           nTicsRadial,                      &
                                           nTicsAngular,                     &
                                           deltaRadial,                      &
                                           deltaAngular,                     &
                                           firstTicRadial,                   &
                                           firstTicAngular,                  &
                                           totalGridPoints,                  &
                                                                  gridWeight )

       real,              intent (in)  :: semiAxisMajor
       real,              intent (in)  :: semiAxisMinor
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadiusMajor
       real,              intent (in)  :: gaussianRadiusMinor
       real,              intent (in)  :: gaussianCenterMajor
       real,              intent (in)  :: gaussianCenterMinor
       integer,           intent (in)  :: nTicsRadial
       integer,           intent (in)  :: nTicsAngular
       real,              intent (in)  :: deltaRadial
       real,              intent (in)  :: deltaAngular
       real,              intent (in)  :: firstTicRadial
       real,              intent (in)  :: firstTicAngular
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam3DGridWeightRadial
  end interface

  interface
     subroutine ed_beam3DGridWeightRecBeam (crossSectionFunctionType,         &
                                            gaussianExponent,                 &
                                            gaussianRadiusMajor,              &
                                            gaussianRadiusMinor,              &
                                            gaussianCenterMajor,              &
                                            gaussianCenterMinor,              &
                                            nTicsSemiAxisMajor,               &
                                            nTicsSemiAxisMinor,               &
                                            deltaSemiAxisMajor,               &
                                            deltaSemiAxisMinor,               &
                                            totalGridPoints,                  &
                                                                   gridWeight )

       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadiusMajor
       real,              intent (in)  :: gaussianRadiusMinor
       real,              intent (in)  :: gaussianCenterMajor
       real,              intent (in)  :: gaussianCenterMinor
       integer,           intent (in)  :: nTicsSemiAxisMajor
       integer,           intent (in)  :: nTicsSemiAxisMinor
       real,              intent (in)  :: deltaSemiAxisMajor
       real,              intent (in)  :: deltaSemiAxisMinor
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam3DGridWeightRecBeam
  end interface

  interface
     subroutine ed_beam3DGridWeightSquare (semiAxisMajor,                    &
                                           semiAxisMinor,                    &
                                           crossSectionFunctionType,         &
                                           gaussianExponent,                 &
                                           gaussianRadiusMajor,              &
                                           gaussianRadiusMinor,              &
                                           gaussianCenterMajor,              &
                                           gaussianCenterMinor,              &
                                           nTicsSemiAxisMajor,               &
                                           nTicsSemiAxisMinor,               &
                                           delta,                            &
                                           totalGridPoints,                  &
                                                                  gridWeight )

       real,              intent (in)  :: semiAxisMajor
       real,              intent (in)  :: semiAxisMinor
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadiusMajor
       real,              intent (in)  :: gaussianRadiusMinor
       real,              intent (in)  :: gaussianCenterMajor
       real,              intent (in)  :: gaussianCenterMinor
       integer,           intent (in)  :: nTicsSemiAxisMajor
       integer,           intent (in)  :: nTicsSemiAxisMinor
       real,              intent (in)  :: delta
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam3DGridWeightSquare
  end interface

  interface
     subroutine ed_beam3DGridWeightStatistical (semiAxisMajor,                    &
                                                semiAxisMinor,                    &
                                                crossSectionFunctionType,         &
                                                gaussianExponent,                 &
                                                gaussianRadiusMajor,              &
                                                gaussianRadiusMinor,              &
                                                gaussianCenterMajor,              &
                                                gaussianCenterMinor,              &
                                                seed,                             &
                                                totalGridPoints,                  &
                                                                       gridWeight )

       real,              intent (in)  :: semiAxisMajor
       real,              intent (in)  :: semiAxisMinor
       character (len=*), intent (in)  :: crossSectionFunctionType
       real,              intent (in)  :: gaussianExponent
       real,              intent (in)  :: gaussianRadiusMajor
       real,              intent (in)  :: gaussianRadiusMinor
       real,              intent (in)  :: gaussianCenterMajor
       real,              intent (in)  :: gaussianCenterMinor
       integer,           intent (in)  :: seed
       integer,           intent (in)  :: totalGridPoints
       real,              intent (out) :: gridWeight
     end subroutine ed_beam3DGridWeightStatistical
  end interface

  interface
     subroutine ed_beamsCheck1DRec () 
     end subroutine ed_beamsCheck1DRec
  end interface

  interface
     subroutine ed_beamsCheck2DCyl3D () 
     end subroutine ed_beamsCheck2DCyl3D
  end interface

  interface
     subroutine ed_beamsCheck2DRec () 
     end subroutine ed_beamsCheck2DRec
  end interface

  interface
     subroutine ed_beamsCheck3DRec () 
     end subroutine ed_beamsCheck3DRec
  end interface

  interface
     subroutine ed_beamsCheck () 
     end subroutine ed_beamsCheck
  end interface

  interface
     subroutine ed_beamsInfo1DRec () 
     end subroutine ed_beamsInfo1DRec
  end interface

  interface
     subroutine ed_beamsInfo2DRec () 
     end subroutine ed_beamsInfo2DRec
  end interface

  interface
     subroutine ed_beamsInfo3DRec () 
     end subroutine ed_beamsInfo3DRec
  end interface

  interface
     subroutine ed_beamsInfo () 
     end subroutine ed_beamsInfo
  end interface

  interface
     subroutine ed_blockData1DRec (blockID,              &
                                   iminBlock, imaxBlock, &
                                   iminData , imaxData,  &
                                   iminDerv , imaxDerv,  &
                                   deltaInvI,            &
                                   blockData             )

       integer, intent (in) :: blockID
       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: iminData , imaxData
       integer, intent (in) :: iminDerv , imaxDerv
       real,    intent (in) :: deltaInvI
       real,    intent (in) :: blockData (:,:)
     end subroutine ed_blockData1DRec
  end interface

  interface
     subroutine ed_blockData2DRec (blockID,              &
                                   iminBlock, imaxBlock, &
                                   jminBlock, jmaxBlock, &
                                   iminData , imaxData,  &
                                   jminData , jmaxData,  &
                                   iminDerv , imaxDerv,  &
                                   jminDerv , jmaxDerv,  &
                                   deltaI   , deltaJ,    &
                                   deltaInvI, deltaInvJ, &
                                   blockData             )

       integer, intent (in) :: blockID
       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: jminBlock, jmaxBlock
       integer, intent (in) :: iminData , imaxData
       integer, intent (in) :: jminData , jmaxData
       integer, intent (in) :: iminDerv , imaxDerv
       integer, intent (in) :: jminDerv , jmaxDerv
       real,    intent (in) :: deltaI   , deltaJ
       real,    intent (in) :: deltaInvI, deltaInvJ
       real,    intent (in) :: blockData (:,:,:)
     end subroutine ed_blockData2DRec
  end interface

  interface
     subroutine ed_blockData3DRec (iminBlock, imaxBlock,            &
                                   jminBlock, jmaxBlock,            &
                                   kminBlock, kmaxBlock,            &
                                   iminData , imaxData,             &
                                   jminData , jmaxData,             &
                                   kminData , kmaxData,             &
                                   iminDerv , imaxDerv,             &
                                   jminDerv , jmaxDerv,             &
                                   kminDerv , kmaxDerv,             &
                                   deltaI   , deltaJ   , deltaK,    &
                                   deltaInvI, deltaInvJ, deltaInvK, &
                                   blockData                        )

       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: jminBlock, jmaxBlock
       integer, intent (in) :: kminBlock, kmaxBlock
       integer, intent (in) :: iminData , imaxData
       integer, intent (in) :: jminData , jmaxData
       integer, intent (in) :: kminData , kmaxData
       integer, intent (in) :: iminDerv , imaxDerv
       integer, intent (in) :: jminDerv , jmaxDerv
       integer, intent (in) :: kminDerv , kmaxDerv
       real,    intent (in) :: deltaI   , deltaJ   , deltaK
       real,    intent (in) :: deltaInvI, deltaInvJ, deltaInvK
       real,    intent (in) :: blockData (:,:,:,:)
     end subroutine ed_blockData3DRec
  end interface

  interface
     subroutine ed_computeBeamPower (timeStep, timeSimulation, pulseID, beamPower)
       real,    intent (in)  :: timeStep
       real,    intent (in)  :: timeSimulation
       integer, intent (in)  :: pulseID
       real,    intent (out) :: beamPower
     end subroutine ed_computeBeamPower
  end interface

  interface
     subroutine ed_checkReuseDepo (reuseDepo,laserIsOn)
       implicit none
       logical, intent(OUT) :: reuseDepo
       logical, intent(in)  :: laserIsOn
     end subroutine ed_checkReuseDepo
  end interface

  interface
     real function ed_CoulombFactor (Z,e,k,T,Ne)
       real, intent (in) :: Z
       real, intent (in) :: e
       real, intent (in) :: k
       real, intent (in) :: T
       real, intent (in) :: Ne
     end function ed_CoulombFactor
  end interface

  interface
     subroutine ed_createRays1DRec (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ed_createRays1DRec
  end interface

  interface
     subroutine ed_createRays2DCyl3D (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ed_createRays2DCyl3D
  end interface

  interface
     subroutine ed_createRays2DRec (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ed_createRays2DRec
  end interface

  interface
     subroutine ed_createRays3DRec (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ed_createRays3DRec
  end interface

  interface
     subroutine ed_createRays (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ed_createRays
  end interface
  
  interface
     subroutine ed_createRayTags (passSplitDriver)
       integer, optional, intent (in) :: passSplitDriver
     end subroutine ed_createRayTags
  end interface

  interface
     subroutine ed_gridEllipseRadial (nGridPoints,               &
                                                   nTicsRadial,  &
                                                   nTicsAngular, &
                                                   deltaRadial,  &
                                                   deltaAngular  )

       integer, intent (inout) :: nGridPoints
       integer, intent (out)   :: nTicsRadial
       integer, intent (out)   :: nTicsAngular
       real,    intent (out)   :: deltaRadial
       real,    intent (out)   :: deltaAngular
     end subroutine ed_gridEllipseRadial
  end interface

  interface
     subroutine ed_gridEllipseSquare (aspectRatio,                     &
                                                   nGridPoints,        &
                                                   nTicsSemiaxisMajor, &
                                                   nTicsSemiaxisMinor, &
                                                   deltaNormalized     )

       real,    intent (in)    :: aspectRatio
       integer, intent (inout) :: nGridPoints
       integer, intent (out)   :: nTicsSemiaxisMajor
       integer, intent (out)   :: nTicsSemiaxisMinor
       real,    intent (out)   :: deltaNormalized
     end subroutine ed_gridEllipseSquare
  end interface

  interface
     subroutine ed_gridLineRegular (lineLength, nGridPoints,    delta)
       real,    intent (in)  :: lineLength
       integer, intent (in)  :: nGridPoints
       real,    intent (out) :: delta
     end subroutine ed_gridLineRegular
  end interface

  interface
     subroutine ed_gridRectangle (aspectRatio,                       &
                                               nGridPoints,          &
                                               nTicsSemiaxisMajor,   &
                                               nTicsSemiaxisMinor,   &
                                               deltaNormalizedMajor, &
                                               deltaNormalizedMinor  )

       real,    intent (in)    :: aspectRatio
       integer, intent (inout) :: nGridPoints
       integer, intent (out)   :: nTicsSemiaxisMajor
       integer, intent (out)   :: nTicsSemiaxisMinor
       real,    intent (out)   :: deltaNormalizedMajor
       real,    intent (out)   :: deltaNormalizedMinor
     end subroutine ed_gridRectangle
  end interface

  interface
     subroutine ed_initializeBeams () 
     end subroutine ed_initializeBeams
  end interface

  interface
     subroutine ed_initializeRays () 
     end subroutine ed_initializeRays
  end interface

  interface
     real function ed_inverseBremsstrahlungRate (Z,e,Me,k,T,Ne,Nc,lnLambda)
       real, intent (in) :: Z
       real, intent (in) :: e
       real, intent (in) :: Me
       real, intent (in) :: k
       real, intent (in) :: T
       real, intent (in) :: Ne
       real, intent (in) :: Nc
       real, intent (in) :: lnLambda
     end function ed_inverseBremsstrahlungRate
  end interface

  interface
     subroutine ed_laserIOCheckWrite (passSplitDriver)
       integer, optional, intent (in) :: passSplitDriver
     end subroutine ed_laserIOCheckWrite
  end interface

  interface
     subroutine ed_laserIOInit ()
     end subroutine ed_laserIOInit
  end interface

  interface
    function ed_maxConfinement1DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_maxConfinement1DRec (1:nc)
    end function ed_maxConfinement1DRec
  end interface

  interface
    function ed_maxConfinement2DCyl3D (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_maxConfinement2DCyl3D (1:nc)
    end function ed_maxConfinement2DCyl3D
  end interface

  interface
    function ed_maxConfinement2DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_maxConfinement2DRec (1:nc)
    end function ed_maxConfinement2DRec
  end interface

  interface
    function ed_maxConfinement3DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_maxConfinement3DRec (1:nc)
    end function ed_maxConfinement3DRec
  end interface

  interface
    function ed_minConfinement1DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_minConfinement1DRec (1:nc)
    end function ed_minConfinement1DRec
  end interface

  interface
    function ed_minConfinement2DCyl3D (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_minConfinement2DCyl3D (1:nc)
    end function ed_minConfinement2DCyl3D
  end interface

  interface
    function ed_minConfinement2DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_minConfinement2DRec (1:nc)
    end function ed_minConfinement2DRec
  end interface

  interface
    function ed_minConfinement3DRec (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ed_minConfinement3DRec (1:nc)
    end function ed_minConfinement3DRec
  end interface

  interface
     subroutine ed_printBeamsData () 
     end subroutine ed_printBeamsData
  end interface

  interface
     subroutine ed_printBlockVariable (blockID, variable, fileUnit)
       integer, intent (in) :: blockID
       integer, intent (in) :: variable
       integer, intent (in) :: fileUnit
     end subroutine ed_printBlockVariable
  end interface

  interface
     subroutine ed_printEnergyStamp (timeStep, timeSimulation) 
       real, intent (in) :: timeStep
       real, intent (in) :: timeSimulation
     end subroutine ed_printEnergyStamp
  end interface

  interface
     subroutine ed_printMainData () 
     end subroutine ed_printMainData
  end interface

  interface
     subroutine ed_printMatrix (fileUnit,                   &
                                title,                      &
                                rowMinMatrix, rowMaxMatrix, &
                                colMinMatrix, colMaxMatrix, &
                                rowMinPrint , rowMaxPrint,  &
                                colMinPrint , colMaxPrint,  &
                                matrix                      )

      integer,           intent (in) :: fileUnit
      character (len=*), intent (in) :: title
      integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
      integer,           intent (in) :: colMinMatrix, colMaxMatrix
      integer,           intent (in) :: rowMinPrint , rowMaxPrint
      integer,           intent (in) :: colMinPrint , colMaxPrint
      real,              intent (in) :: matrix (rowMinMatrix : rowMaxMatrix , colMinMatrix : colMaxMatrix)
     end subroutine ed_printMatrix
  end interface

  interface
     subroutine ed_printPulsesData () 
     end subroutine ed_printPulsesData
  end interface

  interface
     subroutine ed_printRaysData (processorID) 
       integer, intent (in) :: processorID
     end subroutine ed_printRaysData
  end interface

  interface
     subroutine ed_raysBlockIDInfo () 
     end subroutine ed_raysBlockIDInfo
  end interface

  interface
    function ed_raytraceODEfunction1DRec (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: ed_raytraceODEfunction1DRec (1:size (y))
    end function ed_raytraceODEfunction1DRec
  end interface

  interface
    function ed_raytraceODEfunction2DCyl3D (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: ed_raytraceODEfunction2DCyl3D (1:size (y))
    end function ed_raytraceODEfunction2DCyl3D
  end interface

  interface
    function ed_raytraceODEfunction2DRec (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: ed_raytraceODEfunction2DRec (1:size (y))
    end function ed_raytraceODEfunction2DRec
  end interface

  interface
    function ed_raytraceODEfunction3DRec (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: ed_raytraceODEfunction3DRec (1:size (y))
    end function ed_raytraceODEfunction3DRec
  end interface

  interface
     subroutine ed_saveRays () 
     end subroutine ed_saveRays
  end interface

  interface
     subroutine ed_setupBeams () 
     end subroutine ed_setupBeams
  end interface

  interface
     subroutine ed_setupPulses () 
     end subroutine ed_setupPulses
  end interface

  interface
     subroutine ed_setupRays () 
     end subroutine ed_setupRays
  end interface

  interface
     real function ed_time2FacesParabolicPath1D (pos, vel, acc, minFace, maxFace, noFaceTime)
       real, intent (in) :: pos, vel, acc
       real, intent (in) :: minFace, maxFace
       real, intent (in) :: noFaceTime
     end function ed_time2FacesParabolicPath1D
  end interface

  interface
     subroutine ed_traceBlockRays1DRec (timeStep,                          &
                                        rayFirst,  rayLast,                &
                                        iminBlock, imaxBlock,              &
                                        xminBlock, xmaxBlock,              &
                                        deltaX,                            &
                                        deltaInvX,                         &
                                        blockReflectMinX,                  &
                                        blockReflectMaxX,                  &
                                                           cellEnergyDepot ) 

       real,    intent (in)    :: timeStep
       integer, intent (in)    :: rayFirst,  rayLast   
       integer, intent (in)    :: iminBlock, imaxBlock
       real,    intent (in)    :: xminBlock, xmaxBlock
       real,    intent (in)    :: deltaX
       real,    intent (in)    :: deltaInvX
       logical, intent (in)    :: blockReflectMinX
       logical, intent (in)    :: blockReflectMaxX
       real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock)
     end subroutine ed_traceBlockRays1DRec
  end interface

  interface
     subroutine ed_traceBlockRays2DCyl3D (timeStep,                          &
                                          rayFirst,  rayLast,                &
                                          iminBlock, imaxBlock,              &
                                          jminBlock, jmaxBlock,              &
                                          xminBlock, xmaxBlock,              &
                                          zminBlock, zmaxBlock,              &
                                          deltaX, deltaZ,                    &
                                          deltaInvX, deltaInvZ,              &
                                          blockReflectMinX,                  &
                                          blockReflectMaxX,                  &
                                          blockReflectMinZ,                  &
                                          blockReflectMaxZ,                  &
                                                            wedgeEnergyDepot ) 

       real,    intent (in)    :: timeStep
       integer, intent (in)    :: rayFirst,  rayLast   
       integer, intent (in)    :: iminBlock, imaxBlock
       integer, intent (in)    :: jminBlock, jmaxBlock
       real,    intent (in)    :: xminBlock, xmaxBlock
       real,    intent (in)    :: zminBlock, zmaxBlock
       real,    intent (in)    :: deltaX, deltaZ
       real,    intent (in)    :: deltaInvX, deltaInvZ
       logical, intent (in)    :: blockReflectMinX
       logical, intent (in)    :: blockReflectMaxX
       logical, intent (in)    :: blockReflectMinZ
       logical, intent (in)    :: blockReflectMaxZ
       real,    intent (inout) :: wedgeEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)
     end subroutine ed_traceBlockRays2DCyl3D
  end interface

  interface
     subroutine ed_traceBlockRays2DRec (timeStep,                          &
                                        rayFirst,  rayLast,                &
                                        iminBlock, imaxBlock,              &
                                        jminBlock, jmaxBlock,              &
                                        xminBlock, xmaxBlock,              &
                                        yminBlock, ymaxBlock,              &
                                        deltaX, deltaY,                    &
                                        deltaInvX, deltaInvY,              &
                                        blockReflectMinX,                  &
                                        blockReflectMaxX,                  &
                                        blockReflectMinY,                  &
                                        blockReflectMaxY,                  &
                                                           cellEnergyDepot ) 

       real,    intent (in)    :: timeStep
       integer, intent (in)    :: rayFirst,  rayLast   
       integer, intent (in)    :: iminBlock, imaxBlock
       integer, intent (in)    :: jminBlock, jmaxBlock
       real,    intent (in)    :: xminBlock, xmaxBlock
       real,    intent (in)    :: yminBlock, ymaxBlock
       real,    intent (in)    :: deltaX, deltaY
       real,    intent (in)    :: deltaInvX, deltaInvY
       logical, intent (in)    :: blockReflectMinX
       logical, intent (in)    :: blockReflectMaxX
       logical, intent (in)    :: blockReflectMinY
       logical, intent (in)    :: blockReflectMaxY
       real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)
     end subroutine ed_traceBlockRays2DRec
  end interface

  interface
     subroutine ed_traceBlockRays3DRec (timeStep,                          &
                                        rayFirst,  rayLast,                &
                                        iminBlock, imaxBlock,              &
                                        jminBlock, jmaxBlock,              &
                                        kminBlock, kmaxBlock,              &
                                        xminBlock, xmaxBlock,              &
                                        yminBlock, ymaxBlock,              &
                                        zminBlock, zmaxBlock,              &
                                        deltaX, deltaY, deltaZ,            &
                                        deltaInvX, deltaInvY, deltaInvZ,   &
                                        blockReflectMinX,                  &
                                        blockReflectMaxX,                  &
                                        blockReflectMinY,                  &
                                        blockReflectMaxY,                  &
                                        blockReflectMinZ,                  &
                                        blockReflectMaxZ,                  &
                                                          cellEnergyDepot, &
                                                        cellIntensityDepot ) 

       real,    intent (in)    :: timeStep
       integer, intent (in)    :: rayFirst,  rayLast   
       integer, intent (in)    :: iminBlock, imaxBlock
       integer, intent (in)    :: jminBlock, jmaxBlock
       integer, intent (in)    :: kminBlock, kmaxBlock
       real,    intent (in)    :: xminBlock, xmaxBlock
       real,    intent (in)    :: yminBlock, ymaxBlock
       real,    intent (in)    :: zminBlock, zmaxBlock
       real,    intent (in)    :: deltaX, deltaY, deltaZ
       real,    intent (in)    :: deltaInvX, deltaInvY, deltaInvZ
       logical, intent (in)    :: blockReflectMinX
       logical, intent (in)    :: blockReflectMaxX
       logical, intent (in)    :: blockReflectMinY
       logical, intent (in)    :: blockReflectMaxY
       logical, intent (in)    :: blockReflectMinZ
       logical, intent (in)    :: blockReflectMaxZ
       real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock)
       real,    intent (inout),OPTIONAL &
                               :: cellIntensityDepot (iminBlock:      ,jminBlock:         ,kminBlock:)
     end subroutine ed_traceBlockRays3DRec
  end interface

  interface
     subroutine ed_traceRays1DRec (timeStep)
       real, intent (in) :: timeStep
     end subroutine ed_traceRays1DRec
  end interface

  interface
     subroutine ed_traceRays2DCyl3D (timeStep)
       real, intent (in) :: timeStep
     end subroutine ed_traceRays2DCyl3D
  end interface

  interface
     subroutine ed_traceRays2DRec (timeStep)
       real, intent (in) :: timeStep
     end subroutine ed_traceRays2DRec
  end interface

  interface
     subroutine ed_traceRays3DRec (timeStep)
       real, intent (in) :: timeStep
     end subroutine ed_traceRays3DRec
  end interface

  interface
     subroutine ed_traceRays (timeStep)
       real, intent (in) :: timeStep
     end subroutine ed_traceRays
  end interface

  interface
     subroutine ed_updatePlasma (blockCount,blockList,scaleFact)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,OPTIONAL, intent (in) :: scaleFact
     end subroutine ed_updatePlasma
  end interface

  interface
     subroutine ed_updateRays (doMove) 
       logical, intent (in) :: doMove
     end subroutine ed_updateRays
  end interface

end Module ed_interface
