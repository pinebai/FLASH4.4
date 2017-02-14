!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      module io

        integer, save :: iu_log = 6

        public :: output_dir, amr_log_file
#ifdef LIBRARY
        character (len=80) :: output_dir
#else
        character (len=80) :: output_dir = OUTPUT_DIR
#endif
        character (len=80) :: amr_log_file

      end module io
