!! An interface for the MA28 routine.  The internal workings are ignored as a kernel

Module la_interface
  interface 
     subroutine ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag)
       integer, intent(IN)    ::  n,nz,licn,lirn 
       real, intent(IN)       ::  u
       integer, intent(INOUT) ::  irn(lirn),icn(licn),iw(n,8) 
       real, intent(INOUT)    ::  a(licn),w(n)
       integer, intent(OUT)   ::  ikeep(n,5),iflag 
     end subroutine ma28ad
  end interface
end Module la_interface
