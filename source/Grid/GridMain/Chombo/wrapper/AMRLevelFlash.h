#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELFLASH_H_
#define _AMRLEVELFLASH_H_

// DEV CD: Added "QuadCFInterp.H" and "PiecewiseLinearFillPatch.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "QuadCFInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "UsingNamespace.H"

#include "flash_ctypes.h"
typedef std::map <int,Vector<DataIndex> > map_data_index_t;

/// AMR Godunov
/**
 */
class AMRLevelFlash: public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelFlash(int a_level);

  /// Destructor
  /**
   */
  virtual ~AMRLevelFlash();

  /// Define the parameters the object needs
  /**
   */
  void defineParams(const flash_amr_info_t& flashAMRInfo, const mesh_info_t& meshInfo, bool restart);

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep
  /**
   */
  virtual Real advance();

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags);

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags);

  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  // call this instead of initalGrid if this was loaded from a restart
  void restartGrid();

  /// Initialize data
  /**
   */
  virtual void initialData();

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;

  ///
  bool allDefined() const;


  // DEV CD: Things to do after regrid
  virtual void postRegrid(int a_base_level);

  // DEV CD: New public methods added for FLASH
  void zeroFluxData();

  void putFluxData(const void *pFluxes,
		   const Real dt,
		   const int levelBoxID,
		   const int axis);

  void reflux();

  void fillGuardCells();

  void coarseAverage();

  int getNumBoxes(const int procID,
		  const int gridStruct) const;

  const DisjointBoxLayout& getDisjointBoxLayout (const int gridStruct) const;

  DataIndex getDataIndex(const int levelBoxID,
			 const int gridStruct) const;

  void * getDataPtr(const int levelBoxID,
		    const int gridStruct);

  box_info_t getBoxInfo(const int levelBoxID,
			const int gridStruct) const;

  void makeBoxLookupTable();
  

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelFlash* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelFlash* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  // Grid spacing
  Real m_dx;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterp;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  // Number of conserved states
  int m_numStates;

  // Names of conserved states
  Vector<string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // Grid layout for this level
  DisjointBoxLayout m_grids;

  // True if all the parameters for this object are defined
  bool m_paramsDefined;


  // DEV CD: New data members
  map_data_index_t m_boxLookupTable;
  flash_amr_info_t m_flashAMRInfo;
  mesh_info_t m_meshInfo;
  LevelData<FArrayBox> m_scratchCtrFab;

  bool m_useFluxCorrection;
  bool m_useQuadCFInterp;
  bool m_scaleFineFluxes;
  bool m_restart;
private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelFlash& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelFlash(const AMRLevelFlash& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
