!!****f* source/IO/IO_readUserArray
!!
!!  NAME
!!    IO_readUserArray
!!
!!  SYNOPSIS
!!    call IO_readUserArray()
!!
!!
!!  DESCRIPTION
!!
!!    This is the supplied interface for users to read in additional
!!    quantities to the checkpoint or plotfile.  This routine should be used
!!    for reading in various types of arrays.  If the array is a global quantity
!!    only the master processor needs to read in the data.  If it is a quantity 
!!    which is different on all processors then each processor must read in its 
!!    own section of the array. (For a serial IO implementation each processor would
!!    need to send its data to the master.)  The specific implementation is left up
!!    to the user.  
!!
!!    In each case the user should make a call to either 
!!    io_h5read_generic_int_arr (hdf5) or io_ncmpi_read_generic_iarr (pnetcdf)  or
!!    io_h5read_generic_real_arr (hdf5) or io_ncmpi_read_generic_darr (pnetcdf)
!!    depending on the io implementation.
!!  
!!  ARGUMENTS
!!    
!!   
!!
!!  NOTES 
!!
!!    This routine should NOT
!!    be used to read in grid scope data or to read in single scalar
!!    values.  To read in user defined grid scope variables the user should
!!    use the keyword 'GRIDVAR' to declare a grid scope variable in the Config
!!    files.  Then set the runtime parameters plot_grid_var_1, plot_grid_var_2,   
!!    to the name of the grid var to include them in the checkpoint files and
!!    plotfiles.
!!
!!    To read in single scalar quantities the use the IO_setScalar routine to
!!    add a scalar to the scalar output list.
!!
!!
!!  SEE ALSO
!!
!!    io_h5read_generic_int_arr
!!    io_h5read_generic_real_arr
!!    IO_setScalar
!!    
!!    For the pnetcdf implementation see
!!    io_ncmpi_read_generic_iarr
!!    io_ncmpi_read_generic_darr
!!
!!***


subroutine IO_readUserArray ()
  
  implicit none

  

end subroutine IO_readUserArray
