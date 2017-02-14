Module ut_sortInterface
  implicit none  
  interface
     subroutine ut_sortOnProcs(count, props, attrib, &
          numProcs, storage, workspace, &
          perProc, ifNonZeroProc, nonZeroProcsCount)

       implicit none
       
       integer, intent(IN) :: count, props, attrib, numProcs
       real,dimension(props,count),intent(INOUT) ::storage,workspace
       integer,dimension(numProcs),intent(OUT) :: perProc, ifNonZeroProc
       integer, intent(OUT) :: nonZeroProcsCount
       
     end subroutine ut_sortOnProcs
  end interface
end Module ut_sortInterface
