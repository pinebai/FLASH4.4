!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_raytraceODEfunctionData
!!
!! NAME
!!
!!  ed_raytraceODEfunctionData
!!
!! SYNOPSIS
!!
!!  use ed_raytraceODEfunctionData
!!  
!! DESCRIPTION
!!
!!  Data module serving as data transmitter to the ray tracing ODE functions
!!  ------------------------------------------------------------------------
!!   
!!   acc(X,Y,Z)         : acceleration component in (x,y,z)-direction
!!   accFactor(X,Y,Z)   : acceleration multiplying factor for (x,y,z)-direction
!!   cellEdgeInv(X,Y,Z) : inverse of the cell edge length in (x,y,z)-direction
!!   cellZbar           : the cell Zbar value
!!   i,j,k              : the cell locator indices inside a block
!!   Nele (1:4)         : number of electron density (1) and (x,y,z) derivatives (2,3,4)
!!   rayCritDens        : the ray's critical Density
!!   saveComputations   : logical keyword to save computations
!!   wedgeEdgeInv(X,Z)  : inverse of the wedge edge length in (x,z)-direction (3D in 2D cylindrical only)
!!   wedgeSlope         : the wedge slope value (3D in 2D cylindrical only)
!!   wedgeZbar          : the wedge Zbar value (3D in 2D cylindrical only)
!!   (x,y,z)01          : rescaled [0,1] ray (x,y,z) coordinates
!!   (x,y,z)maxCell     : the highest (x,y,z)-coordinate of the cell
!!   (x,y,z)minCell     : the lowest (x,y,z)-coordinate of the cell
!!   (x,z)maxWedge      : the highest (x,z)-coordinate of the wedge (3D in 2D cylindrical only)
!!   (x,z)minWedge      : the lowest (x,z)-coordinate of the wedge (3D in 2D cylindrical only)
!!  
!! NOTES
!!
!!  All of the variables in here need to have the 'threadprivate' omp attribute!
!!
!!***

Module ed_raytraceODEfunctionData

  implicit none

  logical, save :: saveComputations

  integer, save :: i,j,k

  real,    save :: accX, accY, accZ
  real,    save :: accFactorX, accFactorY, accFactorZ
  real,    save :: cellEdgeInvX, cellEdgeInvY, cellEdgeInvZ
  real,    save :: cellZbar
  real,    save :: rayCritDens
  real,    save :: wedgeEdgeInvX, wedgeEdgeInvZ
  real,    save :: wedgeSlope
  real,    save :: wedgeZbar
  real,    save :: x01, y01, z01
  real,    save :: xmaxCell, ymaxCell, zmaxCell
  real,    save :: xminCell, yminCell, zminCell
  real,    save :: xmaxWedge, zmaxWedge
  real,    save :: xminWedge, zminWedge

  real,    save :: Nele (1:4)

  !$omp threadprivate (i,j,k,                                    &
  !$omp                accX, accY, accZ,                         &
  !$omp                accFactorX, accFactorY, accFactorZ,       &
  !$omp                cellEdgeInvX, cellEdgeInvY, cellEdgeInvZ, &
  !$omp                cellZbar,                                 &
  !$omp                rayCritDens,                              &
  !$omp                wedgeEdgeInvX, wedgeEdgeInvZ,             &
  !$omp                wedgeSlope,                               &
  !$omp                wedgeZbar,                                &
  !$omp                x01, y01, z01,                            &
  !$omp                xmaxCell, ymaxCell, zmaxCell,             &
  !$omp                xminCell, yminCell, zminCell,             &
  !$omp                xmaxWedge, zmaxWedge,                     &
  !$omp                xminWedge, zminWedge,                     &
  !$omp                Nele                                      )

end Module ed_raytraceODEfunctionData
