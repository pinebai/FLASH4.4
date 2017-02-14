#ifndef CHOMBO_ADAPTIVE_GRID_H
#define CHOMBO_ADAPTIVE_GRID_H
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <string>
#include "LevelData.H"
#include "FArrayBox.H"
#include "SPMD.H"
#include "DisjointBoxLayout.H"
#include "AMRLevelFlashFactory.h"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "Tuple.H"

typedef struct {
  int level;
  int levelBoxID;
} box_level_info_t;

class Chombo_Adaptive_Grid {

private:

  // The following variables are used by the member functions:
  // IsInitialRefinementDone, BuildInitialGrid, RefineInitialGrid and
  // FinalizeInitialGrid.  These functions along with FLASH's
  // Simulation_initBlock and Grid_markRefineDerefeine collectively behave
  // like AMR::initialGrid.  The variable names are kept the same as
  // AMR::initialGrid, but the variable scope is sometimes different.
  Vector<Vector<Box> > m_old_grids, m_new_grids;
  Vector<IntVectSet> m_old_tags, m_tags;
  int m_top_level, m_max_level, m_finest_level;

  //AMRLevelFlashFactory m_amrFlashFactory;
  Vector<box_level_info_t> m_boxLevelInfo;

  // A value of 0 means no limit to block size (see BRMeshRefine in user guide)
  int m_max_grid_size;
  int m_max_base_grid_size;

  int m_blockFactor;
  Vector<AMRLevel*> m_amrlevels;
  Vector<int> m_ref_ratios;
  BRMeshRefine m_mesh_refine;
  Vector<Vector<Box> > m_amr_grids;
  int m_verbosity;
  bool m_restrictBeforeGhostExchange;

public:
  Chombo_Adaptive_Grid();
  ~Chombo_Adaptive_Grid();

  // Initializes the Chombo Grid
  void Define(const flash_amr_info_t& flashAMRInfo,
	      const mesh_info_t& meshInfo, bool restart);

  // Copied from AMR class
  void makeBaseLevelMesh(Vector<Box>& a_grids) const;

  bool IsInitialRefinementDone() const;
  void BuildInitialGrid();
  void RefineInitialGrid();
  void FinalizeInitialGrid();

  void Regrid(const int a_base_level);

  // Get the blkIDs belonging to this MPI process
  std::vector<int> GetBlockIDs(const int blkType, const int level) const;

  // Get a vector of cell coordinates
  // -- depends on GetBoxInfo()
  std::vector<double> GetCellCoords(const int blkID, const int axis, const int edge, const int guardcell) const;

  // Return an interoperable data structure containing
  // all of the box information required by FLASH
  box_info_t GetBoxInfo(const int blkID, const int gridDataStruct) const;

  void UpdateBoxLevelInfo();

  // Get a pointer to the requested blkID
  void * GetDataPtr(const int blkID, const int gridDataStruct);

  // Fill guard cells as and when FLASH requires
  void FillGuardCells();

  // Average level data after each FLASH update
  void AverageLevelData();

  void ZeroFluxData();

  void PutFluxData(const void *pFluxes, const Real dt, const int blkID, const int axis);

  void Reflux();

  void PostTimeStep();

  void ReadCheckpoint(const char *filename, double &simTime, double &dt, HDF5HeaderData &scalars, HDF5HeaderData &runparms);
  void WriteCheckpoint(const char filename[], double simTime, double dt, const HDF5HeaderData &scalars, const HDF5HeaderData &runparms) const;

  void PrintBoxes(string msg, const Vector<Vector<Box> > &vvb) const;
};

#endif
