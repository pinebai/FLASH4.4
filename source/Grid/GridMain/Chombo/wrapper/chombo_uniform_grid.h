#ifndef CHOMBO_UNIFORM_GRID_H
#define CHOMBO_UNIFORM_GRID_H
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
#include "AMR.H"
#include "constants.h"
#include "flash_ctypes.h"

class Chombo_Uniform_Grid {

private:
  LevelData<FArrayBox> m_ctrFab, m_scratchCtrFab;
  DisjointBoxLayout m_ugBoxLayout;
  ProblemDomain m_computationalDomain;
  Real m_deltas[SpaceDim];
  IntVect m_ghostVect;
  flash_ug_info_t m_flashUGInfo;
  mesh_info_t m_meshInfo;
  int m_verbosity;

public:
  Chombo_Uniform_Grid();
  ~Chombo_Uniform_Grid();
  void Define(const flash_ug_info_t &flashUGInfo,
	      const mesh_info_t &meshInfo, bool restart);
  void FillGuardCells();
  
  void ReadCheckpoint(
    const char *filename,
    double &simTime,
    double &dt,
    HDF5HeaderData &scalars,
    HDF5HeaderData &runparms);
  void WriteCheckpoint(
    const char filename[],
    double simTime,
    double dt,
    const HDF5HeaderData &scalars,
    const HDF5HeaderData &runparms) const;
  
  void* GetDataPtr(const int blkID,
		    const int gridDataStruct);
  std::vector<int> GetBlockIDs(const int blkType,
			       const int level) const;
  std::vector<double> GetCellCoords(const int blkID,
				    const int axis,
				    const int edge,
				    const int guardcell) const;
  box_info_t GetBoxInfo(const int blkID,
			const int gridStruct) const;
  bool IsInitialRefinementDone() const;
  void BuildInitialGrid();
  void RefineInitialGrid();
  void Regrid(const int baseLevel);
  void AverageLevelData();
  void ZeroFluxData();
  void PutFluxData(const void *pFluxes,
		   const Real dt,
		   const int blkID,
		   const int axis);
  void FinalizeInitialGrid();
  void Reflux();
  void PostTimeStep();
  void PrintMeshNames() const;
};

#endif
