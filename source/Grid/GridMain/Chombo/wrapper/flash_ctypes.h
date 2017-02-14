#ifndef FLASH_CTYPES
#define FLASH_CTYPES
#include <map>
#include <vector>
#include <string>
#include "constants.h"
#include "flash_bool.h"

// Interoperable types must be "plain old data structures (PODS)"
// It is essential that they are laid out in memory like C structs.

// flash_ug_info_t is interoperable with Fortran
// TODO: Replace MAX_GRID_DATA_STRUCT_TMP with MAX_GRID_DATA_STRUCT.
// At the moment MAX_GRID_DATA_STRUCT is 5 and so I don't want
// to change 5 to 9 and possibly introduce other problems at this stage.
typedef struct {
  double lowDomain[MDIM];
  double highDomain[MDIM];
  int procGrid[MDIM];
  int baseDomainSize[MDIM];
  int guardCells[MDIM];
  int domainBC[MDIM];
  int meshTypes[MAX_GRID_DATA_STRUCT_TMP];
  int meshNumVars[MAX_GRID_DATA_STRUCT_TMP];
  int verbosity;
} flash_ug_info_t;

// flash_amr_info_t is interoperable with Fortran
typedef struct {
  double lowDomain[MDIM];
  double highDomain[MDIM];
  double BRMeshRefineFillRatio;
  int BRMeshRefineBufferSize;
  int BRMeshRefineBlockFactor;
  int maxBlockSize[MDIM];
  int baseDomainSize[MDIM];
  int guardCells[MDIM];
  int domainBC[MDIM];
  int meshTypes[MAX_GRID_DATA_STRUCT_TMP];
  int meshNumVars[MAX_GRID_DATA_STRUCT_TMP];
  int maxRefineLevel;
  int verbosity;
  int quadCFInterp;
  int fluxCorrect;
  int refRatio;
  int restrictBeforeGhostExchange;
  int scaleFineFluxes;
} flash_amr_info_t;

// box_info_t is interoperable with Fortran
typedef struct {
  double deltas[MDIM];
  double bsize[MDIM];
  double coord[MDIM];
  double lowBndBox[MDIM];
  double highBndBox[MDIM];
  int lowLimits[MDIM];
  int highLimits[MDIM];
  int guardcells[MDIM];
  int corner[MDIM];
  int stride[MDIM];
  int lrefine;
  int isNextToLowDomain[MDIM];
  int isNextToHighDomain[MDIM];
} box_info_t;

typedef struct {
    int real_count;
    char (*real_names)[MAX_STRING_LENGTH];
    double *real_vals;
    
    int int_count;
    char (*int_names)[MAX_STRING_LENGTH];
    int *int_vals;
    
    int str_count;
    char (*str_names)[MAX_STRING_LENGTH];
    char (*str_vals)[MAX_STRING_LENGTH];
    
    int log_count;
    char (*log_names)[MAX_STRING_LENGTH];
    int *log_vals;
} named_vals_t;

// mesh_info_t is not interoperable with Fortran (as intended)
typedef std::map <int,std::vector<std::string> > mesh_info_t;

#endif
