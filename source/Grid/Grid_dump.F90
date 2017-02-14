!!****f* source/Grid/Grid_dump
!!
!! NAME
!!  Grid_dump
!!
!! SYNOPSIS
!!
!!  call Grid_dump(integer(IN) :: var(num),
!!                 integer(IN) :: num,
!!                 integer(IN) :: blockID,
!!                 logical(IN) :: gcell)
!!  
!! DESCRIPTION 
!!  
!! Dumps the variables specified in integer array "var" to a file.
!! This routine doesn't require special resources, so it can be done from   
!! anywhere in the code, and is useful for diagnostic purposes
!! This function can only be used on a single block per processor,
!! and mostly works with Uniform Grid. 
!!  
!! ARGUMENTS 
!!
!!  var :: array containing the indices of the variables to be dumped
!!  num :: number of variables being dumped.
!!  blockID :: local number of block to be dumped
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLE
!!  
!!  num = 3  !dumping 3 variables
!!  var(1) = DENS_VAR
!!  var(2) = PRES_VAR
!!  var(3) = TEMP_VAR
!!  blockID = 1  ! local block number
!!  gcell = .false.
!!
!!  call Grid_dump(var, num, blockID, gcell)
!!  
!!  will dump the interior cell values of density, pressure and temperature
!!  for local block number 1.
!!
!!  To explain the use of gcells, consider a global domain with 8x8 points
!!  mapped on 2x2 processors. each processor has blocks of size
!!  4x4.If there are 2 guard cells along each dimension, then the
!!  block size including guardcells is 8x8 and the distribution on
!!  four processors is as shown below ("*"
!!  are the interior points and the )"o" are guard cells.
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!  If gcell is true then dump is of size 16x16 and looks like
!!
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!
!!  and if gcell is false the dump is of size 8x8 and looks like
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!    
!!
!! NOTES
!!  DENS_VAR, PRES_VAR, TEMP_VAR etc are #defined values in Flash.h
!!  indicating the index in the physical data array.
!!  The routine calling Grid_dump will need to include Flash.h 
!!
!!***

subroutine Grid_dump(var,num,blockID,gcell)

implicit none
  integer, intent(IN) :: num, blockID
  integer, dimension(num), intent(IN) :: var
  logical, intent(IN) :: gcell

  return
end subroutine Grid_dump
