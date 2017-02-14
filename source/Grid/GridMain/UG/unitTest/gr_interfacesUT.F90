Module gr_interfacesUT
#include "constants.h"
  interface
     subroutine gr_testBoundary(dataBlk,strt,fin,guard,axis,bcType,face,isFace)
       implicit none
       integer,intent(IN) :: guard,axis,bcType,face
       integer,dimension(MDIM),intent(IN) :: strt,fin
       logical,intent(IN ) :: isFace
       real,dimension(:,:,:),intent(INOUT) :: dataBlk
     end subroutine gr_testBoundary
  end interface
end Module gr_interfacesUT
