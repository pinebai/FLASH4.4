!!****ih* source/flashUtilities/Pipeline/localAPI/pl_interface
!!
!! NAME
!!
!!  pl_interface
!!
!! SYNOPSIS
!!
!!   use pl_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Pipeline unit that defines its
!!  private interfaces.
!!
!!***

Module pl_interface

  interface
     subroutine pl_localPipelineSetup () 
     end subroutine pl_localPipelineSetup
  end interface

  interface
     subroutine pl_localPostRecvMsg (channel)
       integer, intent (in) :: channel
     end subroutine pl_localPostRecvMsg
  end interface

  interface
     subroutine pl_localPostSendMsg (channel)
       integer, intent (in) :: channel
     end subroutine pl_localPostSendMsg
  end interface

  interface
     subroutine pl_localSaveRecvItems (channel, isSaved)
       integer, intent (in)  :: channel
       logical, intent (out) :: isSaved
     end subroutine pl_localSaveRecvItems
  end interface

  interface
     subroutine pl_printGlobalStatusVector (fileUnit)
       integer, intent (in) :: fileUnit
     end subroutine pl_printGlobalStatusVector
  end interface

end Module pl_interface
