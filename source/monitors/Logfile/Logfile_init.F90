!!****f* source/monitors/Logfile/Logfile_init
!!
!! NAME
!!  Logfile_init
!!
!! SYNOPSIS
!!  Logfile_init()
!!               
!!               
!!
!! DESCRIPTION
!!  Performs Logfile initializations.  Initializes the runtime
!!  runtime parameters needed by the Logfile Unit and then calls
!!  Logfile_create
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used by Logfile unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!    log_file [STRING]
!!        Name of log file to create
!!    run_comment [STRING]
!!        Comment for run
!!    run_number [STRING]
!!        Identification number for run
!!    basenm
!!
!!
!!***


subroutine Logfile_init()
implicit none
end subroutine Logfile_init
