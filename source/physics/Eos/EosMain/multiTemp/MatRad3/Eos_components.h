#if 0
This file contains definitions for using Eos with multiple fluid components.
#endif


#define EOSCOMP_BEGIN 1
#define EOSCOMP_NUM_COMPONENTS N_EOS_TEMP

#define EOSCOMP_ION 1
#define EOSCOMP_ELE 2
#define EOSCOMP_RAD 3

#if 0
The following is not counted in EOSCOMP_NUM_COMPONENTS; it is used in some
places to indicate "normal matter", i.e., "not radiation".
#endif
#define EOSCOMP_MATTER 2
