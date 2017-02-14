!!****ih* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picInterface
!!
!! This is the header file for the initialization subunit of the Particles
!! unit. It defines the interfaces with subunit scope.
!! 
!!***

Module pt_picInterface

  interface
     subroutine pt_picUpdatePositions(dt, halfDt)
       real, intent(in)    :: dt
       logical, intent(in) :: halfDt
     end subroutine pt_picUpdatePositions
  end interface


  interface
     subroutine pt_picEfield(ibx, iby, ibz)
       integer, intent(in) :: ibx, iby, ibz
     end subroutine pt_picEfield
  end interface

  interface
     subroutine pt_picAdvanceVel(dt, halfDt)
       real, intent(IN)    :: dt
       logical, intent(IN) :: halfDt
     end subroutine pt_picAdvanceVel
  end interface

  interface
     subroutine pt_picAdvanceB(dt)
       real, intent(IN)    :: dt
     end subroutine pt_picAdvanceB
  end interface
  
  interface
     subroutine pt_picSetCurrents()
     end subroutine pt_picSetCurrents
  end interface

  interface
     subroutine pt_picInit()
     end subroutine pt_picInit
  end interface

  interface
     subroutine pt_picCurl(ix, iy, iz, jx, jy, jz, mul, zer)
       integer, intent(in) :: ix, iy, iz, jx, jy, jz
       real,    intent(in) :: mul
       logical, intent(in) :: zer
     end subroutine pt_picCurl
  end interface

  interface
     subroutine pt_picApplyBoundary()
     end subroutine pt_picApplyBoundary
  end interface

end Module pt_picInterface
