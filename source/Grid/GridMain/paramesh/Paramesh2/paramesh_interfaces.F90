

module paramesh_interfaces



  interface
     subroutine amr_close
     end subroutine amr_close
  end interface
  
  
  interface
     subroutine amr_guardcell(mype,iopt,nlayers,idir,maxNodetype_gcWanted)
       integer, intent(in)  ::  idir
       integer, intent(in)  ::  mype,iopt,nlayers
       integer,OPTIONAL, intent(in) :: maxNodetype_gcWanted
     end subroutine amr_guardcell

  end interface
  
  interface
      subroutine amr_initialize
      end subroutine amr_initialize
   end interface
   
   
   interface
      subroutine amr_refine_derefine
      end subroutine amr_refine_derefine
   end interface

   interface
      subroutine amr_prolong(MyPE, iopt, nlayers)
        integer, intent(in) :: MyPE, iopt, nlayers
      end subroutine amr_prolong
   end interface

   interface
      subroutine amr_restrict(MyPE,iopt,iempty)
        integer, intent(in) :: MyPE, iopt, iempty
      end subroutine amr_restrict
   end interface

 interface
    subroutine amr_flux_conserve(MyPE, nsub, axis)
      integer, intent(IN) :: MyPE, nsub, axis
    end subroutine amr_flux_conserve
 end interface

 interface
    subroutine amr_prolong_unk_fun & 
         &     (recv, isg, ioff, joff, koff, jface, mype, isrc)
      use physicaldata
      use tree
      integer :: isg, ioff, joff, koff, jface, mype
      real    :: recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      integer,intent(IN),OPTIONAL :: isrc
    end subroutine amr_prolong_unk_fun
 end interface

 interface
    subroutine amr_prolong_gen_unk_fun & 
         &       (recv, ia, ib, ja, jb, ka, kb, isg, & 
         &        ioff, joff, koff, mype, isrc)
#include "Flash.h"
      real, intent(IN), &
           dimension(NUNK_VARS,GRID_ILO_GC:GRID_IHI_GC,&
           GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC)::recv
      integer, intent(IN) :: isg, mype, ia, ib, ja, jb, ka, kb,ioff,joff,koff
      integer,intent(IN),OPTIONAL :: isrc
    end subroutine amr_prolong_gen_unk_fun
 end interface

 end module paramesh_interfaces

