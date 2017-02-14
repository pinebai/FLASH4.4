!!****f* source/physics/Cosmology/Cosmology_redshiftHydro
!!
!!
!! NAME
!!
!!  Cosmology_redshiftHydro
!!
!! SYNOPSIS
!!
!!  Cosmology_redshiftHydro( integer(IN) :: blockCount, 
!!                           integer(IN) :: blockList(blockCount)) 
!!
!! DESCRIPTION
!!
!!  Description:  Applies the redshift operator to hydrodynamical quantities.
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList -   array holding local IDs of blocks on which to advance
!!
!!
!!***


subroutine Cosmology_redshiftHydro ( blockCount, blockList)
  
  implicit none

  integer, INTENT(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList

end subroutine Cosmology_redshiftHydro

