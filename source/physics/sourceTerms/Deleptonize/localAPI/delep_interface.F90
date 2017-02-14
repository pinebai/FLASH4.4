!!****ih* source/physics/sourceTerms/Deleptonize/localAPI/delep_interface
!!
!!   Deleptonize's private interfaces
!!***

Module delep_interface

  interface
     subroutine delep_detectBounce (blockCount,blockList,dt,time)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt,time
     end subroutine delep_detectBounce
  end interface

end Module delep_interface
