!!****h* source/monitors/Logfile/Logfile_interface
!!
!! This is the header file for the Logfile module
!! that defines its public interfaces.
!!***
Module Logfile_interface

  implicit none

#include "Flash.h"
#include "constants.h"

  interface
    subroutine Logfile_break ( char)
      character(len=1), intent(in)       :: char
    end subroutine Logfile_break
  end interface

  interface
    subroutine Logfile_close(logUnitLocal)
      logical,optional,intent(IN) :: logUnitLocal
    end subroutine Logfile_close
  end interface

  interface
    subroutine Logfile_create ()
    end subroutine Logfile_create
  end interface

  interface
    subroutine Logfile_finalize()
    end subroutine Logfile_finalize
  end interface

  interface
    subroutine Logfile_getDateTimeStr(dateTimeStr)
      implicit none
      character(len=28), intent(OUT)        :: dateTimeStr
    end subroutine Logfile_getDateTimeStr
  end interface

  interface
    subroutine Logfile_init()
     end subroutine Logfile_init
  end interface

  interface
     subroutine Logfile_open(logUnit,logUnitLocal)
       integer,intent(OUT) :: logUnit
       logical,intent(IN) :: logUnitLocal
     end subroutine Logfile_open
  end interface
  
  interface Logfile_stamp
    subroutine Logfile_stampInt(intVal, tag, attrib)
      integer, intent(in)                   ::  intVal  
      character(len=*), intent(in)           :: tag
      character(len=*), intent(in), OPTIONAL :: attrib
    end subroutine Logfile_stampInt

    subroutine Logfile_stampIntArray( intArr, len, tag, attrib)
      character(len=*),intent(in)           :: tag
      character(len=*), intent(in),OPTIONAL :: attrib
      integer,intent(in)                    :: len
      integer, dimension(len), intent(in)   :: intArr
    end subroutine Logfile_stampIntArray

    subroutine Logfile_stampReal( realVal, tag, attrib)
      real, intent(in)                       :: realVal
      character(len=*), intent(in)           :: tag
      character(len=*), intent(in), OPTIONAL  :: attrib
    end subroutine Logfile_stampReal

    subroutine Logfile_stampRealArray( realArr, len, tag, attrib)
      character(len=*),intent(in)            :: tag
      character(len=*), intent(in), OPTIONAL :: attrib
      integer,intent(in)                     :: len
      real, dimension(len),intent(in)        :: realArr
    end subroutine Logfile_stampRealArray

    subroutine Logfile_stampStr( string, tag, attrib)
      character(len=*), intent(in)           :: string
      character(len=*), intent(in), OPTIONAL :: tag, attrib
    end subroutine Logfile_stampStr

    subroutine Logfile_stampStrArray( strArr, len, tag, attrib)
      character(len=*), intent(in)             :: tag
      character(len=*), intent(in), OPTIONAL   :: attrib
      integer,intent(in)                       :: len
      character(len=*), dimension(len),intent(in) :: strArr
    end subroutine Logfile_stampStrArray

    subroutine Logfile_stampStrPair( strArr, len, dimen, tag, attrib)
      character(len=*),intent(in)           :: tag
      character(len=*),intent(in), OPTIONAL :: attrib
      integer,intent(in)                    :: len, dimen
      character(len=*),dimension(len,dimen),intent(in) :: strArr
    end subroutine Logfile_stampStrPair
  end interface

  interface
    subroutine Logfile_stampMessage( string,force)
      character(len=*), intent(in)           :: string
      logical, intent(in), optional          :: force
    end subroutine Logfile_stampMessage
  end interface

  interface
     subroutine Logfile_stampVarMask(unkVarMask, willCallEos, tag, maskTag)
       logical, intent(in) :: unkVarMask(:)
       logical, intent(in) :: willCallEos
       character(len=*),intent(in) :: tag, maskTag
     end subroutine Logfile_stampVarMask
  end interface

  interface
    subroutine Logfile_writeSummary(strArr, length, dim, strLen, numHeaders,reduced,separateFiles)
      integer, intent(in)                      :: length, dim, strLen, numHeaders
      character(len=MAX_STRING_LENGTH), dimension(length,dim), intent(in)  :: strArr
      logical, optional, intent(IN)                        :: reduced
      logical, optional, intent(IN)                        :: separateFiles
    end subroutine Logfile_writeSummary
  end interface

end Module Logfile_interface
