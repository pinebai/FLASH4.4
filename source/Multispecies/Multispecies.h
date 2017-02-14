#if 0
This is the header file for the Multispecies database.  For all routines
several properties are used. They are described below and defined with
an integer at the bottom of the file.
    description             property
    -------------------------------------------------------------------------------------
    numTotal                A             Total number of protons and neutrons in nucleus
    numPositive             Z             Atomic number, number of protons in nucleus;
                                                 upper bound for Z* (ionization level)
    zMin                    MS_ZMIN       Lower bound for Z*
    numNeutral              N             Number of neutrons
    numNegative             E             Number of electrons
    bindingEnergy           EB            Binding energy
    adiabatic index         GAMMA         Ratio of heat capacities: Cp / Cv
    opacity lowTemp limit   MS_OPLOWTEMP  Temperature in K below which the low temperature
                                          opacity calculation kicks in

    EOS type                MS_EOSTYPE    Gamma-based, table, or something else?

#endif

#define A (1)
#define Z (2)
#define N (3)
#define E (4)
#define EB (5)
#define GAMMA 6
#define MS_ZMIN 7
#define MS_EOSTYPE 107
#define MS_EOSSUBTYPE 108
#define MS_EOSZFREEFILE 208
#define MS_EOSENERFILE 209
#define MS_EOSPRESFILE 210
#define MS_EOSGROUPNAME 211
#define MS_EOSIONFILE 408
#define MS_EOSELEFILE 409
#define MS_NUMELEMS   300
#define MS_ZELEMS     301
#define MS_AELEMS     302
#define MS_FRACTIONS  303
#define MS_OPLOWTEMP  304

#define UNDEFINED_REAL -999.
#define UNDEFINED_INT -999
#define MS_UNDEFINED_STRING '-none-'
#define NO_MASK -1

#define MS_STRINGLEN 128
