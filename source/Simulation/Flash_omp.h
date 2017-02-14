#if(0)
 This header file has some definitions that are useful for parts of the FLASH code
 that uses OpenMP directives. It is currently not necessary that all source files
 with OpenMP directives include this file. However, those that make use of the macros
 defined here need to include it.
#endif


#if(0)
! If _OPENMP is defined AND indicates a version of 200505 or earlier, then we assume
! the compiler processes omp directives but would fail on the collapse(2) clause in them.
! We shall have COLLAPSE(2) expand to nothing in that case. That means that OpenMP 'do'
! directives will be applied to the outer loop, which should be good enough at least in
! the cases where it is really worth using OpenMP, i.e., in 3D setups.
! We name the macro COLLAPSE instead of anything else just in case some compiler uses
! a proprocessor that does not expand macros in Fortran comments and/or directives:
! in that case, the unmodified 'COLLAPSE(2)' will be left unchanged by the preprocessor
! to be handled as part of the OpenMP directive.
#endif

#ifdef _OPENMP
#if _OPENMP > 200505
#define COLLAPSE collapse
#else
#define COLLAPSE(x)
#endif
#endif
