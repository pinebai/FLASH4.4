!!****ih* source/diagnostics/ProtonImaging/localAPI/pi_interface
!!
!! NAME
!!
!!  pi_interface
!!
!! SYNOPSIS
!!
!!   use pi_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Proton Imaging unit that defines its
!!  private interfaces.
!!
!!***

Module pi_interface

  interface
     subroutine pi_beamsCheck2DCyl3D () 
     end subroutine pi_beamsCheck2DCyl3D
  end interface

  interface
     subroutine pi_beamsCheck3DRec () 
     end subroutine pi_beamsCheck3DRec
  end interface

  interface
     subroutine pi_beamsCheck () 
     end subroutine pi_beamsCheck
  end interface

  interface
     subroutine pi_beamsInfo3DRec () 
     end subroutine pi_beamsInfo3DRec
  end interface

  interface
     subroutine pi_beamsInfo () 
     end subroutine pi_beamsInfo
  end interface

  interface
     subroutine pi_blockData3DRec (iminBlock, imaxBlock,            &
                                   jminBlock, jmaxBlock,            &
                                   kminBlock, kmaxBlock,            &
                                   deltaInvX, deltaInvY, deltaInvZ, &
                                   blockData                        )

       integer, intent (in) :: iminBlock, imaxBlock
       integer, intent (in) :: jminBlock, jmaxBlock
       integer, intent (in) :: kminBlock, kmaxBlock
       real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
       real,    intent (in) :: blockData (:,:,:,:)
     end subroutine pi_blockData3DRec
  end interface

  interface
     subroutine pi_calculateCellData (coordMin, coordMax, indexMin, indexMax, delta, coords)
       real,    intent (in)  :: coordMin, coordMax
       integer, intent (in)  :: indexMin, indexMax
       real,    intent (out) :: delta
       real,    intent (out) :: coords (indexMin:indexMax+1)
     end subroutine pi_calculateCellData
  end interface

  interface
     subroutine pi_capsuleGrainIndices2xyz (i,j,k,S,x,y,z)
       integer, intent (in)  :: i,j,k
       real,    intent (in)  :: S
       real,    intent (out) :: x,y,z
     end subroutine pi_capsuleGrainIndices2xyz
  end interface

  interface
     subroutine pi_capsuleNextGrainIndices (L,i,j,k,valid)
       integer, intent (in)    :: L
       integer, intent (inout) :: i,j,k
       logical, intent (out)   :: valid
     end subroutine pi_capsuleNextGrainIndices
  end interface

  interface
     integer function pi_capsuleTotalGrainCount (L)
       integer, intent (in) :: L
     end function pi_capsuleTotalGrainCount
  end interface

  interface
     subroutine pi_closeDetectorFiles ()
     end subroutine pi_closeDetectorFiles
  end interface

  interface
     subroutine pi_closeDiskProtonFile (whichOne)
       character (len=3), intent (in) :: whichOne
     end subroutine pi_closeDiskProtonFile
  end interface

  interface
     subroutine pi_closeMonitorFile ()
     end subroutine pi_closeMonitorFile
  end interface

  interface
     subroutine pi_createProtons3DRec (blockCount, blockList, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeSimulation
     end subroutine pi_createProtons3DRec
  end interface

  interface
     subroutine pi_createProtons (blockCount, blockList, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeSimulation
     end subroutine pi_createProtons
  end interface

  interface
     subroutine pi_createProtonTags () 
     end subroutine pi_createProtonTags
  end interface

  interface
     subroutine pi_detectorsCheck3DRec () 
     end subroutine pi_detectorsCheck3DRec
  end interface

  interface
     subroutine pi_detectorsCheck () 
     end subroutine pi_detectorsCheck
  end interface

  interface
     subroutine pi_detectorsInfo3DRec () 
     end subroutine pi_detectorsInfo3DRec
  end interface

  interface
     subroutine pi_detectorsInfo () 
     end subroutine pi_detectorsInfo
  end interface

  interface
     subroutine pi_flushDiskProtons2Disk () 
     end subroutine pi_flushDiskProtons2Disk
  end interface

  interface
     subroutine pi_flushScreenProtons2Disk () 
     end subroutine pi_flushScreenProtons2Disk
  end interface

  interface
     subroutine pi_IOdetectorScreens () 
     end subroutine pi_IOdetectorScreens
  end interface

  interface
     subroutine pi_IOinit () 
     end subroutine pi_IOinit
  end interface

  interface
     subroutine pi_IOprotonsCapsule2Domain () 
     end subroutine pi_IOprotonsCapsule2Domain
  end interface

  interface
     function pi_maxConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: pi_maxConfinement3DRec (1:nc)
     end function pi_maxConfinement3DRec
  end interface

  interface
     function pi_minConfinement3DRec (nc,y)
       integer, intent (in) :: nc
       real,    intent (in) :: y (:)
       real                 :: pi_minConfinement3DRec (1:nc)
     end function pi_minConfinement3DRec
  end interface

  interface
     subroutine pi_openDetectorFiles (timeSimulation) 
       real, intent (in) :: timeSimulation
     end subroutine pi_openDetectorFiles
  end interface

  interface
     subroutine pi_openDiskProtonFile (whichOne, saveOldRecords)
       character (len=3), intent (in) :: whichOne
       logical, optional, intent (in) :: saveOldRecords
     end subroutine pi_openDiskProtonFile
  end interface

  interface
     subroutine pi_openMonitorFile ()
     end subroutine pi_openMonitorFile
  end interface

  interface
     real function pi_parabolicPathLength3D (vx,vy,vz,ax,ay,az,T)
       real, intent (in) :: vx,vy,vz
       real, intent (in) :: ax,ay,az
       real, intent (in) :: T
     end function pi_parabolicPathLength3D
  end interface

  interface
     subroutine pi_printBeamsData () 
     end subroutine pi_printBeamsData
  end interface

  interface
     subroutine pi_printBlockVariable (blockID, variable, fileUnit)
       integer, intent (in) :: blockID
       integer, intent (in) :: variable
       integer, intent (in) :: fileUnit
     end subroutine pi_printBlockVariable
  end interface

  interface
     subroutine pi_printDetectorsData () 
     end subroutine pi_printDetectorsData
  end interface

  interface
     subroutine pi_printMainData () 
     end subroutine pi_printMainData
  end interface

  interface
     subroutine pi_printMatrix (fileUnit,                   &
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
     end subroutine pi_printMatrix
  end interface

  interface
     subroutine pi_printProtonsData (fileLabel, processorID) 
       character (len=*), intent (in) :: fileLabel
       integer,           intent (in) :: processorID
     end subroutine pi_printProtonsData
  end interface

  interface
     subroutine pi_protonsBlockIDInfo () 
     end subroutine pi_protonsBlockIDInfo
  end interface

  interface
    subroutine pi_readDiskProtons (blockCount, blockList, moreOnDisk)
      integer, intent (in)    :: blockCount
      integer, intent (in)    :: blockList (1:blockCount)
      logical, intent (inout) :: moreOnDisk
    end subroutine pi_readDiskProtons
  end interface

   interface
    subroutine pi_readDiskProtons3DRec (blockCount, blockList, moreOnDisk)
      integer, intent (in)    :: blockCount
      integer, intent (in)    :: blockList (1:blockCount)
      logical, intent (inout) :: moreOnDisk
    end subroutine pi_readDiskProtons3DRec
  end interface

 interface
     subroutine pi_recordProtonOnScreen (px, py, pz,                &
                                         vx, vy, vz,                &
                                         Jv, Kx, Ky, Kz,            &
                                         detector,                  &
                                                         onScreen,  &
                                                         sx, sy, sz )
       real,    intent (in)            :: px, py, pz
       real,    intent (in)            :: vx, vy, vz
       real,    intent (in)            :: Jv, Kx, Ky, Kz
       integer, intent (in)            :: detector
       logical, intent (out), optional :: onScreen
       real,    intent (out), optional :: sx, sy, sz
     end subroutine pi_recordProtonOnScreen
  end interface

  interface
     subroutine pi_setupBeams () 
     end subroutine pi_setupBeams
  end interface

  interface
     subroutine pi_setupDetectorFileNames (timeSimulation)
       real, intent (in) :: timeSimulation
     end subroutine pi_setupDetectorFileNames
  end interface

  interface
     subroutine pi_setupDetectors () 
     end subroutine pi_setupDetectors
  end interface

  interface
     subroutine pi_setupDiskProtons () 
     end subroutine pi_setupDiskProtons
  end interface

  interface
     subroutine pi_setupProtons () 
     end subroutine pi_setupProtons
  end interface

  interface
     subroutine pi_setupScreenProtons () 
     end subroutine pi_setupScreenProtons
  end interface

  interface
     subroutine pi_statisticalFinalize ()
     end subroutine pi_statisticalFinalize
  end interface

  interface
     subroutine pi_statisticalInit ()
     end subroutine pi_statisticalInit
  end interface

  interface
     subroutine pi_statisticalSetSeed (seed)
       integer, intent (in) :: seed
     end subroutine pi_statisticalSetSeed
  end interface

  interface
     subroutine pi_statisticalXYcircle (radius,             &
                                        arraySize,          &
                                                   nCircle, &
                                                   xCircle, &
                                                   yCircle  )

       real,    intent (in)  :: radius
       integer, intent (in)  :: arraySize
       integer, intent (out) :: nCircle
       real,    intent (out) :: xCircle (1:arraySize)
       real,    intent (out) :: yCircle (1:arraySize)
     end subroutine pi_statisticalXYcircle
  end interface

  interface
     subroutine pi_statisticalXYZsphere (radius,             &
                                         arraySize,          &
                                                    nSphere, &
                                                    xSphere, &
                                                    ySphere, &
                                                    zSphere  )

       real,    intent (in)  :: radius
       integer, intent (in)  :: arraySize
       integer, intent (out) :: nSphere
       real,    intent (out) :: xSphere (1:arraySize)
       real,    intent (out) :: ySphere (1:arraySize)
       real,    intent (out) :: zSphere (1:arraySize)
     end subroutine pi_statisticalXYZsphere
  end interface

  interface
     real function pi_time2FacesParabolicPath1D (pos, vel, acc, minFace, maxFace, noFaceTime)
       real, intent (in) :: pos, vel, acc
       real, intent (in) :: minFace, maxFace
       real, intent (in) :: noFaceTime
     end function pi_time2FacesParabolicPath1D
  end interface

  interface
     subroutine pi_traceBlockProtons3DRec (protonFirst,  protonLast,          &
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
                                           blockReflectMaxZ                   ) 

       integer, intent (in)    :: protonFirst, protonLast   
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
     end subroutine pi_traceBlockProtons3DRec
  end interface

  interface
     function pi_traceODEfunction3DRec (t,y)
       real, intent (in) :: t
       real, intent (in) :: y (:)
       real :: pi_traceODEfunction3DRec (1:size (y))
     end function pi_traceODEfunction3DRec
  end interface

  interface
     subroutine pi_traceProtons3DRec () 
     end subroutine pi_traceProtons3DRec
  end interface

  interface
     subroutine pi_traceProtons () 
     end subroutine pi_traceProtons
  end interface

  interface
     subroutine pi_transportBeamProtons (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine pi_transportBeamProtons
  end interface

  interface
     subroutine pi_transportDiskProtons (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine pi_transportDiskProtons
  end interface

  interface
     subroutine pi_updateProtons (doMove) 
       logical, intent (in) :: doMove
     end subroutine pi_updateProtons
  end interface

  interface
     subroutine pi_write2DetectorFile (detector, recordCount, bucketCount, bucket)
       integer, intent (in) :: detector
       integer, intent (in) :: recordCount
       integer, intent (in) :: bucketCount
       real,    intent (in) :: bucket (1:recordCount,1:bucketCount)
     end subroutine pi_write2DetectorFile
  end interface

  interface
     subroutine pi_xferDiskProtonsNew2OldFile (rewindOldFile) 
       logical, intent (in) :: rewindOldFile
     end subroutine pi_xferDiskProtonsNew2OldFile
  end interface

end Module pi_interface
