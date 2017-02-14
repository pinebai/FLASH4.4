#if 0
This file contains definitions for using Eos with multiple fluid components.
THIS variant is used in setups that do not actually have 3 componenets.
#endif


#define EOSCOMP_BEGIN 1
#define EOSCOMP_NUM_COMPONENTS N_EOS_TEMP

#define EOSCOMP_ION 1
#define EOSCOMP_ELE 1
#define EOSCOMP_RAD -1

#if 0
The following is not counted in EOSCOMP_NUM_COMPONENTS; it is used in some
places to indicate "normal matter", i.e., "not radiation".
#endif
#define EOSCOMP_MATTER 0
