!!****f* source/Grid/Grid_pfftInit
!!
!! NAME 
!!
!!   Grid_pfftInit
!!
!! SYNOPSIS
!!
!!   Grid_pfftInit(integer(IN)  :: ndim,
!!                 logical(IN)  :: needMap
!!                 integer(IN)  :: globalLen(MDIM),
!!                 integer(OUT) :: localArraySize(MDIM),
!!                 integer(IN),optional :: transformType(MDIM),
!!                 integer(IN),optional :: baseDatType(0:MDIM),
!!                 integer(IN),optional :: jProcs,
!!                 integer(IN),optional :: kProcs,
!!                 integer(IN),optional :: refinementLevel,
!!                 real(IN), optional   :: region_bndBox)
!!
!! DESCRIPTION 
!!  
!!  This is the initialization routine for using Pfft. If needMap
!!  is true,  routine creates a schedule of data transfers that would
!!  map AMR or UG grid to PFFT. This schedule is remembered internally
!!  by Pfft until the Grid_pfftFinalize routine is called. If needMap
!!  is false, it is assumed the data distribution is already compatible
!!  with Pfft requirements: That is along IAXIS, the entire row is 
!!  within one processor.
!!
!!  The routine also calls gr_pfftInitMetaData, which is responsible for
!!  calculating the trignometric tables, creating communicators necessary
!!  for distributed transposes in Pfft and allocates workspace
!!
!! ARGUMENTS
!!
!!  ndim          - dimensionality of the problem
!!  needMap       - should be true if the default shape is not compatible with
!!                  requirements of the input Pfft array. Only if
!!                  this argument is true, is the map determined.
!!  globalLen     -  the globalsize of the domain
!!  localArraySize       - after mapping to pfft_grid, the local size for the domain
!!  transformType -  type of transform along each dimension
!!                   if none is specified we assume realtocomplex
!!                   in first dimension and complextocomplex in the rest
!!  baseDatType -    basic data type (rela or complex) along each dimension,
!!                   after the transform for that direction (if any),
!!                   the 0 component specifies data type before IAXIS transform.
!!  jProcs,kProcs - if they are present, they decide the shape of the 
!!                  processor grid for pfft. If they are not present,
!!                  a routine that can automatically determine the shape is 
!!                  called.
!!  refinementLevel - The block refinement level at which we will create the map.
!!  region_bndBox  - If the map is from a subset of the physical domain in AMR to
!!                   Pfft grid, this argument contains the bounding box of the region
!!
!!***
subroutine Grid_pfftInit( ndim, needMap, globalLen, localArraySize, &
     transformType, baseDatType, jProcs, kProcs, refinementLevel, region_bndBox)

#include "constants.h"

  implicit none
  integer, intent(IN) :: ndim
  logical, intent(IN) :: needMap
  integer, dimension(MDIM), intent(IN) :: globalLen
  integer,dimension(MDIM),intent(OUT) ::  localArraySize
  integer,dimension(MDIM),optional,intent(IN) :: transformType
  integer,dimension(0:MDIM),optional,intent(IN) :: baseDatType
  integer,optional,intent(IN) :: jProcs, kProcs
  integer,optional,intent(IN) :: refinementLevel
  real, dimension(LOW:HIGH, MDIM), optional, intent(IN) :: region_bndBox


  localArraySize=1

  return
end subroutine Grid_pfftInit
