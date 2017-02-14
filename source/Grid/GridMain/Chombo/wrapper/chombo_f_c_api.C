#include "chombo_f_c_api.h"
#include <set>

// My rule is that:
// 1) extern "C" functions (i.e. called from Fortran) have unit-based arguments
// 2) chombo grid methods return zero-based indices

// This means the extern "C" functions are responsible for translating
// between unit-based and zero-based indices

using namespace std;

template<class T>
static inline T* talloc(size_t n) {
    return (T*)malloc(n*sizeof(T));
}

extern "C" void c_free(void **p) {
    free(*p);
    *p = 0;
}

extern "C" void ch_define_uniform_grid(const flash_ug_info_t *pFlashUGInfo,
				       const int meshStringLens[],
				       const char meshStrings[],
                       bool restart)
{
  const flash_ug_info_t flashUGInfo = *pFlashUGInfo;
  mesh_info_t meshInfo;

  assert_valid_simulation();
  extract_mesh_strings(flashUGInfo.meshNumVars, flashUGInfo.meshTypes,
		       meshStringLens, meshStrings, meshInfo);

#ifdef CHOMBO_UG
  if(!pChomboGrid) pChomboGrid = new Chombo_Uniform_Grid();
  pChomboGrid->Define(flashUGInfo, meshInfo, restart);
#else
  pChomboGrid = NULL;
#endif
}


extern "C" void ch_define_adaptive_grid(const flash_amr_info_t *pFlashAMRInfo,
					const int meshStringLens[],
					const char meshStrings[],
                    bool restart)
{
  const flash_amr_info_t flashAMRInfo = *pFlashAMRInfo;
  mesh_info_t meshInfo;

  assert_valid_simulation();
  extract_mesh_strings(flashAMRInfo.meshNumVars, flashAMRInfo.meshTypes,
		       meshStringLens, meshStrings, meshInfo);

#ifdef CHOMBO_AMR
  if(!pChomboGrid) pChomboGrid = new Chombo_Adaptive_Grid();
  pChomboGrid->Define(flashAMRInfo, meshInfo, restart);
#else
  pChomboGrid = NULL;
#endif
}


extern "C" void ch_is_initial_refinement_done(int *pIsRefComplete)
{
  bool bIsRefComplete = pChomboGrid->IsInitialRefinementDone();
  *pIsRefComplete = (bIsRefComplete ? FLASH_TRUE : FLASH_FALSE);
}


extern "C" void ch_build_initial_grid()
{
  pChomboGrid->BuildInitialGrid();
}


extern "C" void ch_refine_initial_grid()
{
  pChomboGrid->RefineInitialGrid();
}


extern "C" void ch_finalize_initial_grid()
{
  pChomboGrid->FinalizeInitialGrid();
}


extern "C" void ch_regrid(const int baseLevel)
{
  const int zeroBaseLevel = baseLevel - 1;
  pChomboGrid->Regrid(zeroBaseLevel);
}


extern "C" void ch_restrict_all_levels()
{
  pChomboGrid->AverageLevelData();
}


extern "C" void ch_fill_guardcells()
{
  pChomboGrid->FillGuardCells();
}

// convert fortran string to c by trimming trailing spaces
static string fstrtoc(const char c[], int n) {
    string s(c, n);
    while(n > 1 && c[n-1] == ' ') n--;
    s.resize(n);
    return s;
}
// convert c string to fortran, pad end with trailing spaces
static void cstrtof(const string &cstr, char *fstr, size_t n) {
    size_t m = cstr.copy(fstr, n);
    memset(fstr + m, ' ', n - m);
}
// split string by delimiting character
// ex: "abc,123,,hi" -> ["abc","123","","hi"]
static vector<string> split(const string &s, char delim) {
    vector<string> ss;
    string u;
    for(int i=0; i < s.size(); i++) {
        if(s[i] != delim)
            u += s[i];
        else {
            ss.push_back(u);
            u.clear();
        }
    }
    ss.push_back(u);
    return ss;
}

static void writeNamedVals(HDF5HeaderData *h, const named_vals_t *nv) {
    string s;
    for(int j=0; j < nv->real_count; j++) {
        s = fstrtoc(nv->real_names[j], sizeof(nv->real_names[j]));
        h->m_real[s] = nv->real_vals[j];
    }
    for(int j=0; j < nv->int_count; j++) {
        s = fstrtoc(nv->int_names[j], sizeof(nv->int_names[j]));
        h->m_int[s] = nv->int_vals[j];
    }
    for(int j=0; j < nv->str_count; j++) {
        s = fstrtoc(nv->str_names[j], sizeof(nv->str_names[j]));
        h->m_string[s] = fstrtoc(nv->str_vals[j], sizeof(nv->str_vals[j]));
    }
    string loglist;
    for(int j=0; j < nv->log_count; j++) {
        s = fstrtoc(nv->log_names[j], sizeof(nv->log_names[j]));
        loglist += s;
        if(j < nv->log_count-1) loglist += ';';
        h->m_int[s] = nv->log_vals[j] ? 1 : 0;
    }
    if(nv->log_count > 0)
        h->m_string["FLASH_Logicals"] = loglist;
}

static void readNamedVals(const HDF5HeaderData *h, named_vals_t *nv) {
    // reals
    nv->real_count = h->m_real.size();
    nv->real_names = talloc<char[MAX_STRING_LENGTH]>(nv->real_count);
    nv->real_vals = talloc<double>(nv->real_count);
    int i = 0;
    for(map<string,Real>::const_iterator it=h->m_real.begin(); it!=h->m_real.end(); ++it) {
        cstrtof(it->first, nv->real_names[i], sizeof(nv->real_names[i]));
        nv->real_vals[i] = it->second;
        i++;
    }
    // figure out which ints are really logicals
    vector<string> loglist;
    if(h->m_string.count("FLASH_Logicals") != 0)
        loglist = split(h->m_string.find("FLASH_Logicals")->second, ';');
    set<string> logset(loglist.begin(), loglist.end());
    // ints
    nv->int_count = h->m_int.size() - logset.size();
    nv->int_names = talloc<char[MAX_STRING_LENGTH]>(nv->int_count);
    nv->int_vals = talloc<int>(nv->int_count);
    i = 0;
    for(map<string,int>::const_iterator it=h->m_int.begin(); it!=h->m_int.end(); ++it) {
        if(logset.count(it->first) != 0) continue; // skip logicals encoded as ints
        cstrtof(it->first, nv->int_names[i], sizeof(nv->int_names[i]));
        nv->int_vals[i] = it->second;
        i++;
    }
    // strings
    nv->str_count = h->m_string.size();
    nv->str_names = talloc<char[MAX_STRING_LENGTH]>(nv->str_count);
    nv->str_vals = talloc<char[MAX_STRING_LENGTH]>(nv->str_count);
    i = 0;
    for(map<string,string>::const_iterator it=h->m_string.begin(); it!=h->m_string.end(); ++it) {
        if(it->first == "FLASH_Logicals")
            nv->str_count--;
        else {
            cstrtof(it->first, nv->str_names[i], sizeof(nv->str_names[i]));
            cstrtof(it->second, nv->str_vals[i], sizeof(nv->str_vals[i]));
            i++;
        }
    }
    // logicals
    nv->log_count = logset.size();
    nv->log_names = talloc<char[MAX_STRING_LENGTH]>(nv->log_count);
    nv->log_vals = talloc<int>(nv->log_count);
    i = 0;
    for(map<string,int>::const_iterator it=h->m_int.begin(); it!=h->m_int.end(); ++it) {
        if(logset.count(it->first) == 0) continue; // skip ints
        cstrtof(it->first, nv->log_names[i], sizeof(nv->log_names[i]));
        nv->log_vals[i] = it->second;
        i++;
    }
}

extern "C" void ch_write_checkpoint(
    const char filename[],
	double simTime,
	double dt,
    const named_vals_t *scalars,
    const named_vals_t *runparms
) {    
    HDF5HeaderData hsc, hrp;
    writeNamedVals(&hsc, scalars);
    writeNamedVals(&hrp, runparms);
    pChomboGrid->WriteCheckpoint(filename, simTime, dt, hsc, hrp);
}

extern "C" void ch_read_checkpoint(
    const char filename[],
    double *simTime,
    double *dt,
    named_vals_t *scalars,
    named_vals_t *runparms
) {
#ifdef CHOMBO_UG
   if(!pChomboGrid) pChomboGrid = new Chombo_Uniform_Grid();
#else
   if(!pChomboGrid) pChomboGrid = new Chombo_Adaptive_Grid();
#endif
    HDF5HeaderData hsc, hrp;
    pChomboGrid->ReadCheckpoint(filename, *simTime, *dt, hsc, hrp);
    readNamedVals(&hsc, scalars);
    readNamedVals(&hrp, runparms);
}

extern "C" void ch_get_blk_ptr(const int blkID,
			       const int gridDataStruct,
			       void **ppv)
{
  const int zeroBlkID = blkID - 1;
  *ppv = pChomboGrid->GetDataPtr(zeroBlkID, gridDataStruct);
}


// Chombo does not impose a limit on the number of boxes on a
// single processor, but FLASH does.  We copy the elements
// from a vector into a statically sized array.
extern "C" void ch_get_block_ids(const int blkType,
				 const int level,
				 int blkIDs[],
				 int *pNumBlks)
{
  const int zeroLevel = level - 1;
  vector<int> vBlk;
  int numBlks;
  vBlk = pChomboGrid->GetBlockIDs(blkType, zeroLevel);
  numBlks = vBlk.size();
  assert (numBlks <= MAXBLOCKS);
  for (int i=0; i<numBlks; ++i)
  {
    blkIDs[i] = vBlk[i] + 1; // Convert to unit-based block ID.
  }
  *pNumBlks = numBlks;
}


extern "C" void ch_get_num_blocks(int *pNumBlks)
{
  vector<int> vBlk;
  const int blkType = ALL_BLKS;
  const int level = -1; // ignored by method
  vBlk = pChomboGrid->GetBlockIDs(blkType, level);
  *pNumBlks = vBlk.size();
}


extern "C" void ch_get_cell_coords(const int blkID,
				   const int axis,
				   const int edge,
				   const int size,
				   const int guardcell,
				   double coords[])
{
  const int zeroBlkID = blkID - 1;
  const int zeroAxis = axis - 1;
  vector<double> vCoords;
  int numCells;
  
  assert(zeroAxis >= 0 && zeroAxis < MDIM);
  if (zeroAxis >= SpaceDim)
  {
    assert(size >= 1);
    coords[0] = 0.0;
  }
  else
  {
    vCoords = pChomboGrid->GetCellCoords(zeroBlkID,zeroAxis,edge,guardcell);
    numCells = vCoords.size();
    assert(size >= numCells); // So we stay within bounds of coords array
    for (int i=0; i<numCells; ++i)
    {
      coords[i] = vCoords[i];
    }
  }
}


extern "C" void ch_get_single_cell_coords(const int blkID,
					  const int edge,
					  const int beginCount,
					  int ind[],
					  double coords[])
{
  const int zeroBlkID = blkID - 1;
  vector<double> vCoords;
  int pos, guardcell;

  assert (beginCount == INTERIOR || beginCount == EXTERIOR);
  if (beginCount == EXTERIOR) guardcell = 1;
  if (beginCount == INTERIOR) guardcell = 0;

  for (int i=0; i<MDIM; ++i)
  {
    if (i < SpaceDim)
    {
      vCoords = pChomboGrid->GetCellCoords(zeroBlkID,i,edge,guardcell);
      pos = ind[i] - 1;
      assert (pos >= 0 && pos < vCoords.size()); // Zero based vector access.
      coords[i] = vCoords[pos];
    }
    else
    {
      coords[i] = 0.0;
    }
  }
}

// Box info includes block limits so it is important to
// specify gridDataStruct.
extern "C" void ch_get_box_info(const int blkID,
				const int gridDataStruct,
				box_info_t *pBoxInfo)
{
  const int zeroBlkID = blkID - 1;
  *pBoxInfo = pChomboGrid->GetBoxInfo(zeroBlkID, gridDataStruct);
  unit_base_box_info(pBoxInfo);
}


extern "C" void ch_zero_flux_data()
{
  pChomboGrid->ZeroFluxData();
}


extern "C" void ch_put_flux_data(const void *pFluxes,
				 const double dt,
				 const int blkID,
				 const int axis)
{
  const int zeroBlkID = blkID - 1;
  const int zeroAxis = axis - 1;
  pChomboGrid->PutFluxData(pFluxes, dt, zeroBlkID, zeroAxis);
}


extern "C" void ch_reflux()
{
  pChomboGrid->Reflux();
}


extern "C" void ch_post_time_step()
{
  pChomboGrid->PostTimeStep();
}


extern "C" void ch_finalize()
{
  delete (pChomboGrid);
}


// Helper functions below this point (i.e. not callable by Fortran).

void assert_valid_simulation()
{
  int err = 0;

  if (NDIM != SpaceDim)
  {
    err = Driver_abortFlashC("FLASH dimensionality is different to Chombo "
			     "dimensionality.\nMake sure you are linking "
			     "FLASH against the correct Chombo library");
  }

  if (sizeof(Real) != sizeof(double))
  {
    err = Driver_abortFlashC("The Chombo library has not promoted its custom "
			     "Real type to double precision as is needed by "
			     "FLASH.\nRebuild Chombo with PRECISION=DOUBLE "
                             "and then rebuild FLASH.");
  }

#ifndef CH_MPI
  err = Driver_abortFlashC("FLASH depends on a parallel installation of Chombo."
			   "\nRebuild Chombo with MPI=TRUE and then "
                           "rebuild FLASH.");
#endif

#ifdef CH_USE_MEMORY_TRACKING
  err = Driver_abortFlashC("Chombo memory tracking does not work with FLASH "
			   "applications.\nRebuild Chombo with USE_MT=FALSE "
                           "and then rebuild FLASH.");
#endif

  if (err != 0) exit(EXIT_FAILURE);
}


void extract_mesh_strings(const int meshNumVars[],
			  const int meshTypes[],
			  const int meshStringLens[],
			  const char meshStrings[],
			  mesh_info_t &meshInfo)
{
  const char *p;
  int numVariables, k, strLen, meshDefineVal;

  p = meshStrings;
  k = 0;
  for (int i=0; i<MAX_GRID_DATA_STRUCT_TMP; ++i)
  {
    vector<string> v;
    numVariables = meshNumVars[i];
    meshDefineVal = meshTypes[i];
    for (int j=0; j<numVariables; ++j)
    {
      // The individual strings are not null-terminated
      strLen = meshStringLens[k];
      v.push_back(string(p, strLen));
      ++k;
      p += strLen;
    }
    meshInfo[meshDefineVal] = v;
    v.clear();
  }
}


void unit_base_box_info(box_info_t *pBoxInfo)
{
  for (int i=0; i<MDIM; ++i)
  {
    pBoxInfo->lowLimits[i] += 1;
    pBoxInfo->highLimits[i] += 1;
    pBoxInfo->corner[i] += 1;
    pBoxInfo->lrefine += 1;
  }
}
