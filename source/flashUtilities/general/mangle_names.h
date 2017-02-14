/****** source/IO/mangle_names.h
 *
 * NAME
 *      FTOC
 * SYNPOPSIS
 *      FTOC(C_function_name)(args....)
 * FUNCTION
 *      This will contain a number of macros for mangling
 *      C routine names into something that FORTRAN can link to.
 * SWITCHES
 *      So far, this deals only with 
 *         NOUNDERSCORE
 *      which is what the IBMs need.
 *
 ******/

#ifdef IBM
#define NOUNDERSCORE
#endif

#ifdef NOUNDERSCORE
#define FTOC(x) x
#else
#define FTOC(x) x##_
#endif
