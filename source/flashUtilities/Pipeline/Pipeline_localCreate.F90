!!****f* source/flashUtilities/Pipeline/Pipeline_localCreate
!!
!! NAME
!!  
!!  Pipeline_localCreate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localCreate (integer,           intent (in)           :: itemSize,
!!                             integer,           intent (in)           :: maxItems,
!!                             integer,           intent (in)           :: channelSize,
!!                             character (len=*), intent (in), optional :: logName)
!!
!! DESCRIPTION
!!
!!  Creates a local pipeline section on a processor inside a FLASH application. A (global)
!!  pipeline is defined as a set of processors connected by a closed system of channels
!!  through which items will 'flow' from one processor to another. Feeding of the pipeline
!!  with items is done by a special feeding API, which is called from the unit where the
!!  pipeline is being used for. As a default, the pipeline structure is defined from the
!!  current structure of the grid, but can also be defined in an abstract way, as long as
!!  the channels between each of the processors are set up in a consistent (closed) way.
!!  This means that no open channels can exist, for which sending takes place but no
!!  corresponding receive is ever posted, a situation which would lead to deadlocks in
!!  communication.
!!
!!  An optional 'logName' can be passed, indicating that the user wants to know details
!!  details about the workings of the pipeline for this particular application. The log
!!  file created is specific for each processor, i.e. each processor writes only to its
!!  own log file with name: <basenm>Pipeline<logName>_<procID>.log , where <basenm> is the
!!  base name of the current application and <procID> is the processor ID of the current
!!  processor in 000x format.
!!
!!  The only arguments the user has to provide is: 1) the size (number of elements) of each
!!  item, 2) the maximum number of such items and 3) the size of each pipeline channel, i.e.
!!  the number of items each channel can maximally hold, before sending them to the
!!  connecting processor. Details about the pipeline structure to be set up is deferred to
!!  an internal routine 'pl_localPipelineSetup', thus giving the user the possibility
!!  for setting up his personal pipeline (i.e. not based on the default -> grid structure)
!!  in case he needs it.
!!
!! ARGUMENTS
!!
!!  itemSize    : size of each item (number of elements per item)
!!  maxItems    : maximum number of items
!!  channelSize : channel size (number of items each channel can hold at a time)
!!  logName     : optional log file name
!!
!! NOTES
!!
!!  1) This operation is local. No other processors are involved.
!!
!!  2) Channels are only defined as processor ID's that are not equal to the current processor
!!     ID. If no such channels are defined, the number of channels is zero.
!!
!!  3) Pipelines can only be used sequentially (one after the other). No simultaneous pipelines
!!     can be in operation right now.
!!
!!***

subroutine Pipeline_localCreate (itemSize, maxItems, channelSize, logName)

  implicit none

  character (len=*), intent (in), optional :: logName
  integer,           intent (in)           :: itemSize
  integer,           intent (in)           :: maxItems
  integer,           intent (in)           :: channelSize

  return
end subroutine Pipeline_localCreate
