!!****f* source/multiprocessorTools/Pipeline/Pipeline_globalCheckStructure
!!
!! NAME
!!  
!!  Pipeline_globalCheckStructure
!!
!! SYNOPSIS
!! 
!!  Pipeline_globalCheckStructure ()
!!
!! DESCRIPTION
!!
!!  Checks the structure of the pipeline for closure. This routine is very
!!  important, as open structured pipelines will lead to possible deadlocks
!!  when unconnected channels send items, but there is no matching receiving
!!  channel at the other end.
!!
!!  Closure of a pipeline is defined by two properties: 1) each processor
!!  must have channels connecting only to other processors, not to itself
!!  and 2) all sending channels must have matching receiving channels on other
!!  processors. Point 1) is a local property and is easy to check. The second
!!  point could in principle be checked in two different ways:
!!
!!    a) Have actual non-blocking sends issued for a token item (an integer,
!!       for example) for each channel and post speculating non-blocking
!!       receives for all channels on each processor. Continuously issue
!!       non-blocking probes (mpi_iprobe calls) to check, if all token items
!!       were received properly or if there are any dangling non-blocking
!!       sends and receives waiting for completion.
!!
!!    b) Construct an actual pipeline map from the info about the channel
!!       processors given. In this case there are no actual items sent.
!!
!!  We chose to implement version b), which requires more memory (2x number
!!  of processors logicals) than a), but is sure to deliver the required
!!  info. Version a) needs less memory, but depends on the implemented mpi
!!  structure. In version a) we don't know beforehand how many calls to
!!  mpi_iprobe we need until we have message saturation (i.e. all messages
!!  that could be delivered were actually delivered).
!!
!!  Version b) is implemented the following way: On each processor define
!!  two logical vectors: a sending channel vector and a receiving channel
!!  vector. The sending channel vector is established locally on each
!!  processor and contains .true. in those places where there are connecting
!!  channels. The receiving channel vector is established non-locally using
!!  the all to all mpi command. Closure of the pipeline is achieved by having
!!  all send channel vector entries match those of the corresponding receive
!!  channel vector entries (n = processor):
!!
!!             sendChannel (n)   |   recvChannel (n)  |  pipeline is
!!            -------------------|--------------------|--------------
!!                   true        |       true         |    closed
!!                  false        |      false         |    closed
!!                   true        |      false         |     open
!!                  false        |       true         |     open
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  The routine is global and acts as a barrier to all processors due to the
!!  all to all communication type.
!!
!!***

subroutine Pipeline_globalCheckStructure ()

  implicit none

  return
end subroutine Pipeline_globalCheckStructure
