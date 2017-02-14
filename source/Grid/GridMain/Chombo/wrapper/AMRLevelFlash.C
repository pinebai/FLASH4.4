#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iomanip>

#include "parstream.H"
#include "ParmParse.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "CH_HDF5.H"
#include "AMRIO.H"
#include "AMRLevel.H"

#include "AMRLevelFlash.h"

// A large amount of code in this class is copied directly from
// example/AMRGodunov/srcPolytropic/AMRLevelPolytropicGas.cpp.  
// Some methods have been customized and some new methods have
// been created for FLASH.  New methods are at the bottom of this class.

// Constructor
AMRLevelFlash::AMRLevelFlash(int a_level) {
    m_level = a_level;
    if(s_verbosity >= 3)
        pout() << "AMRLevelFlash default constructor" << endl;
    m_paramsDefined = false;
}

// Destructor
AMRLevelFlash::~AMRLevelFlash() {
    if(s_verbosity >= 3)
        pout() << "AMRLevelFlash destructor" << endl;
    m_paramsDefined = false;
}

void AMRLevelFlash::defineParams(
    const flash_amr_info_t& flashAMRInfo,
    const mesh_info_t& meshInfo,
    bool restart) {
    m_restart = restart;
    m_flashAMRInfo = flashAMRInfo;
    m_meshInfo = meshInfo;

    // Set the physical dimension of the x-axis side of the domain.
    // We can work out the physical size of the y-axis and z-axis by
    // multiplying the x-axis cell spacing by the number of cells along a
    // specific axis (see Chombo VisIt plugin).  The cell spacing dx is fixed.
    m_domainLength = flashAMRInfo.highDomain[0] - flashAMRInfo.lowDomain[0];

    verbosity(flashAMRInfo.verbosity);

    // Set the number of ghost cells
    m_numGhost = flashAMRInfo.guardCells[0];
    for (int i=0; i<SpaceDim; ++i) {
        CH_assert(m_numGhost == flashAMRInfo.guardCells[i]);
    }

    mesh_info_t::const_iterator it = meshInfo.find(CENTER);
    CH_assert (it != meshInfo.end());
    if(!restart) {
        m_stateNames = it->second;
        m_numStates = m_stateNames.size();
    }
    else {
        // perhaps do some validation?
    }
    
    m_useFluxCorrection = (FLASH_TRUE == flashAMRInfo.fluxCorrect);
    m_useQuadCFInterp = (FLASH_TRUE == flashAMRInfo.quadCFInterp);
    m_scaleFineFluxes = (FLASH_TRUE == flashAMRInfo.scaleFineFluxes);

    m_paramsDefined = true;
}

// This instance should never get called - historical
void AMRLevelFlash::define(
    AMRLevel*  a_coarserLevelPtr,
    const Box& a_problemDomain,
    int a_level,
    int a_refRatio
) {
    MayDay::Error("AMRLevelFlash::define:\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelFlash::define(
    AMRLevel* a_coarserLevelPtr,
    const ProblemDomain& a_problemDomain,
    int a_level,
    int a_refRatio
) {
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::define " << a_level << endl;
  }

  // DEV CD: I'm assuming that the defineParams method is called by
  // AMRLevelFlashFactory to set parameters before this define method is called
  CH_assert(m_paramsDefined == true);

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelFlash* amrGodPtr = dynamic_cast<AMRLevelFlash*>(a_coarserLevelPtr);

    if (amrGodPtr != NULL)
    {
      // DEV CD: See m_paramsDefined assertion above
    }
    else
    {
      MayDay::Error("AMRLevelFlash::define: a_coarserLevelPtr is not castable to AMRLevelFlash*");
    }
  }

  // DEV CD: Every Chombo example that I see assumes that there is a
  // fixed cell spacing for each axis.  If we wish to break this
  // assumption and use a non-unit aspect ratio then (in addition to
  // making sure FLASH Hydro works with a non-unit aspect ratio) we need
  // to add extra attributes to the Chombo checkpoint file.
  //
  // The Chombo plugin for VisIt checks for attributes named "prob_lo"
  // and "aspect_ratio" in file avtChomboFileFormat.C.  These attributes
  // allow the physical lower boundary to have a different value than
  // 0.0 and the possibility for cells to have a non-unit aspect ratio.
  // Despite VisIt checking for the existance of these attributes, I
  // cannot find any reference to these attributes in Chombo-3.0.
  // Perhaps VisIt Chombo plugin developers were trying to think 
  // one step ahead of the Chombo developers???  In any case, if we need
  // to create applications without current restrictions then we should
  // add the just described attributes and VisIt will work as expected.
  //
  // Cell spacing = dx (The same for each axis)
  // Domain size = sx, sy (Possibly different for x-axis and y-axis)
  // Number of cells = nx, ny (Possibly different for x-axis and y-axis)
  // Assert that the user specifies a grid that leads to a fixed dx.
  // We use the following expression to avoid division rounding errors:
  // sx / nx = dx, sy / ny = dx, and so (sx.ny) = (sy.nx)
  
  if (SpaceDim >= 2)
  {
    const IntVect vn = a_problemDomain.domainBox().size();
    const RealVect vs
      (D_DECL6(
	 m_flashAMRInfo.highDomain[0] - m_flashAMRInfo.lowDomain[0],
	 m_flashAMRInfo.highDomain[1] - m_flashAMRInfo.lowDomain[1],
	 m_flashAMRInfo.highDomain[2] - m_flashAMRInfo.lowDomain[2],
	 0,0,0));

    for (int i=1; i<SpaceDim; ++i)
    {
      const Real sn1 = vs[0] * vn[i];
      const Real sn2 = vs[i] * vn[0];

      if (sn1 != sn2)
      {
	MayDay::Error("AMRLevelFlash::define: non-unit aspect ratio");
      }
    }
  }


  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().size(0);

  CH_assert(isDefined());

  // DEV CD: m_numGhost, m_numStates, m_stateNames initialized in defineParams
}

// Advance by one timestep
Real AMRLevelFlash::advance()
{
  CH_assert(allDefined());

  // We want to do the advancement in FLASH's Hydro subroutine and not here
  MayDay::Error("AMRLevelFlash::advance: Should not be called by FLASH");

  return 0.0;
}

// Things to do after a timestep
void AMRLevelFlash::postTimeStep()
{
  CH_assert(allDefined());

  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::postTimeStep " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // DEV CD: Reflux call moved into a separate method named reflux
    // so that FLASH can correct fluxes after each split PPM sweep.

    // Average from finer level data
    AMRLevelFlash* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                    amrGodFinerPtr->m_UNew);
  }

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelFlash::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Interval curComp(comp,comp);
      Real integral = computeSum(m_UNew,NULL,nRefFine,m_dx,curComp);

      pout() << "    " << setw(23)
             << setprecision(16)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific)
             << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 && !first)
      {
        pout() << " (" << setw(23)
               << setprecision(16)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << (integral-last_integral)
               << " " << setw(23)
               << setprecision(16)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << (integral-orig_integral)
               << ")";
      }

      pout() << endl;

      if (comp == 0)
      {
        if (first)
        {
          orig_integral = integral;
          first = false;
        }

        last_integral = integral;
      }
    }
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::postTimeStep " << m_level << " finished" << endl;
  }
}

// DEV CD: FLASH updates m_scratchCtrFab data with FLASH_TRUE or
// FLASH_FALSE in gr_markRefineDerefine to indicate whether or not a
// cell is tagged.  The method AMRLevelFlash::tagCells adds all cells
// that are marked with FLASH_TRUE.  It is important that FLASH performs
// the interpolation to undefined ghost cells and a ghost cell exchange
// before tagging cells. This happens in Grid_markRefineDerefine subroutine
// by calling Grid_fillGuardCells (and thus AMRLevelFlash::fillGuardCells)
// before calling gr_markRefineDerefine.

// Create tags for regridding
void AMRLevelFlash::tagCells(IntVectSet& a_tags)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::tagCells " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::tagCellsInit " << m_level << endl;
  }

  // Find the mesh variable index for tagging by matching the string "tagc".
  const int notFound = -1;
  int tag_index = notFound;
  CH_assert (m_meshInfo.find(SCRATCH_CTR) != m_meshInfo.end());
  const Vector<std::string> &v = m_meshInfo[SCRATCH_CTR];

  for (int j=0; j<v.size(); ++j)
  {
    if (v[j] == "tagc")
    {
      tag_index = j;
    }
  }

  // If found then read each cell's "tagc" variable and if it is FLASH_TRUE (i.e. the
  // cell is marked for refinement) then add the tagged cell to the IntVectSet.
  if (notFound == tag_index)
  {
    if (s_verbosity >= 3)
      pout() << "Did not find a mesh variable for tagging" << endl;
  }
  else
  {
    if (s_verbosity >= 3)
      pout() << "The mesh variable for tagging is named " <<
	v[tag_index] << " and is at index " << tag_index << endl;

    // I copied the code for iterating over all *internal* cells from 
    // test/quadCFInterpTest.cpp.
    DataIterator dit = m_scratchCtrFab.dataIterator();
    const DisjointBoxLayout& dbl = m_scratchCtrFab.disjointBoxLayout();

    // Iterate over each box in the FArrayBox layout.
    for(dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& farray = m_scratchCtrFab[dit()];

      // Iterate over the internal cells of each box.
      // Note: farray.box() would return a box including guardcells.
      const Box& internalBox = dbl.get(dit());
      BoxIterator bit(internalBox);
      for(bit.reset(); bit.ok(); ++bit)
      {
	//pout() << "AMRLevelFlash::tagCells " << bit() << " = " <<
	//farray(bit(),tag_index) << endl;
	const IntVect& iv = bit();
	if (FLASH_TRUE == farray(bit(),tag_index))
	{
	  a_tags |= iv;
	}
      }
    }
  }
}

// Create tags at initialization
void AMRLevelFlash::tagCellsInit(IntVectSet& a_tags)
{
  CH_assert(allDefined());

  tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelFlash::regrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::regrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;

    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Save data for later
  DataIterator dit = m_UNew.dataIterator();
  for (; dit.ok(); ++dit)
  {
    m_UOld[dit()].copy(m_UNew[dit()]);
  }

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost * IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelFlash* amrGodCoarserPtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
  }

  // Copy from old state
  m_UOld.copyTo(m_UOld.interval(),
                m_UNew,
                m_UNew.interval());

  m_UOld.define(m_grids,m_numStates,ivGhost);
}

void AMRLevelFlash::restartGrid() {
    levelSetup();
}

// Initialize grids
void AMRLevelFlash::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::initialGrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelFlash::initialData()
{
  CH_assert(allDefined());

  // DEV CD: We initialize data in FLASH's Simulation_initBlock.
  // initialData method should not be called.
  MayDay::Error("AMRLevelFlash::initialData: Should not be called by FLASH");
}

// Things to do after initialization
void AMRLevelFlash::postInitialize()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::postInitialize " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Volume weighted average from finer level data
    AMRLevelFlash* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                    amrGodFinerPtr->m_UNew);
  }
}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelFlash::writeCheckpointHeader(HDF5Handle& a_handle) const {
    CH_assert(allDefined());
    if(s_verbosity >= 3) pout() << "AMRLevelFlash::writeCheckpointHeader" << endl;

    HDF5HeaderData header;
    // Setup the number of components
    header.m_int["num_components"] = m_numStates;
    // Setup the component names
    for(int comp = 0; comp < m_numStates; ++comp) {
        char compStr[30];
        sprintf(compStr, "component_%d", comp);
        header.m_string[compStr] = m_stateNames[comp];
    }
    // Write the header
    header.writeToFile(a_handle);
    if(s_verbosity >= 3) pout() << header << endl;
}

// Write checkpoint data for this level
void AMRLevelFlash::writeCheckpointLevel(HDF5Handle& a_handle) const {
    CH_assert(allDefined());
    if(s_verbosity >= 3) pout() << "AMRLevelFlash::writeCheckpointLevel" << endl;

    { // enter the level group
        char levelStr[30];
        sprintf(levelStr, "level_%d", m_level);
        a_handle.setGroup(levelStr);
    }

    // Setup the level header information
    HDF5HeaderData header;
    header.m_int ["ref_ratio"] = m_ref_ratio;
    header.m_real["dx"] = m_dx;
    header.m_real["dt"] = m_dt;
    header.m_real["time"] = m_time;
    header.m_box ["prob_domain"] = m_problem_domain.domainBox();
    // Setup the periodicity info
    D_TERM(
        header.m_int["is_periodic_0"] = m_problem_domain.isPeriodic(0) ? 1 : 0;,
        header.m_int["is_periodic_1"] = m_problem_domain.isPeriodic(1) ? 1 : 0;,
        header.m_int["is_periodic_2"] = m_problem_domain.isPeriodic(2) ? 1 : 0;
    );
    // Write the header
    header.writeToFile(a_handle);
    if(s_verbosity >= 3) pout() << header << endl;

    // Write the data for this level
    write(a_handle, m_UNew.boxLayout());
    write(a_handle, m_UNew, "data");
}

// Read checkpoint header
void AMRLevelFlash::readCheckpointHeader(HDF5Handle& a_handle) {
    if(s_verbosity >= 3) pout() << "AMRLevelFlash::readCheckpointHeader" << endl;

    // Reader the header
    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if(s_verbosity >= 3)
        pout() << "hdf5 header data:" << endl << header << endl;

    // Get the number of components
    if(header.m_int.find("num_components") == header.m_int.end())
        MayDay::Error("AMRLevelFlash::readCheckpointHeader: checkpoint file does not have num_components");
    m_numStates = header.m_int["num_components"];

    // Get the component names
    m_stateNames.resize(m_numStates);
    for(int comp = 0; comp < m_numStates; ++comp) {
        char compStr[60];
        sprintf(compStr, "component_%d", comp);
        if(!header.m_string.count(compStr))
            MayDay::Error("AMRLevelFlash::readCheckpointHeader: checkpoint file does not have enough component names");
        m_stateNames[comp] = header.m_string[compStr];
    }
}

// Read checkpoint data for this level
// define() and restartGrid() must be called after this
void AMRLevelFlash::readCheckpointLevel(HDF5Handle& a_handle) {
    if(s_verbosity >= 3) pout() << "AMRLevelFlash::readCheckpointLevel" << endl;
    
    { // Set the level group
        char levelStr[30];
        sprintf(levelStr, "level_%d", m_level);
        a_handle.pushGroup(levelStr);
    }

    HDF5HeaderData header;
    header.readFromFile(a_handle);
    if(s_verbosity >= 3) pout() << "hdf5 header data:" << endl << header << endl;

    if(!header.m_int.count("ref_ratio"))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain ref_ratio");
    m_ref_ratio = header.m_int["ref_ratio"];
    if(s_verbosity >= 2) pout() << "read ref_ratio = " << m_ref_ratio << endl;

    if(!header.m_real.count("dx"))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain dx");
    m_dx = header.m_real["dx"];
    if(s_verbosity >= 2) pout() << "read dx = " << m_dx << endl;

    if(!header.m_real.count("dt"))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain dt");
    m_dt = header.m_real["dt"];
    if(s_verbosity >= 2) pout() << "read dt = " << m_dt << endl;

    if(!header.m_real.count("time"))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain time");
    m_time = header.m_real["time"];
    if(s_verbosity >= 2) pout() << "read time = " << m_time << endl;

    if(!header.m_box.count("prob_domain"))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain prob_domain");
    Box domainBox = header.m_box["prob_domain"];

    // Get the periodicity info -- this is more complicated than it really
    // needs to be in order to preserve backward compatibility
    bool isPeriodic[SpaceDim];
    D_TERM(
        isPeriodic[0] = header.m_int.count("is_periodic_0") && header.m_int["is_periodic_0"] == 1;,
        isPeriodic[1] = header.m_int.count("is_periodic_1") && header.m_int["is_periodic_1"] == 1;,
        isPeriodic[2] = header.m_int.count("is_periodic_2") && header.m_int["is_periodic_2"] == 1;
    );
    m_problem_domain = ProblemDomain(domainBox, isPeriodic);

    // Get the grids
    if(0 != read(a_handle, m_level_grids))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain a Vector<Box>");
    
    // Create level domain
    m_grids = loadBalance(m_level_grids);
    if(s_verbosity >= 4) {
        pout() << "read level domain: " << endl;
        LayoutIterator lit = m_grids.layoutIterator();
        for(lit.begin(); lit.ok(); ++lit) {
            const Box& b = m_grids[lit()];
            pout() << lit().intCode() << ": " << b << endl;
        }
        pout() << endl;
    }
    
    // Reshape state with new grids
    m_UNew.define(m_grids, m_numStates);
    if(0 != read<FArrayBox>(a_handle, m_UNew, "data", m_grids))
        MayDay::Error("AMRLevelFlash::readCheckpointLevel: file does not contain state data");
    m_UOld.define(m_grids, m_numStates);
    
    a_handle.popGroup();
}

// Write plotfile header
void AMRLevelFlash::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelFlash::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelFlash::writePlotLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());

  // DEV CD: Pass IntVect::Zero to not print guard cells.
  write(a_handle,m_UNew,"data", IntVect::Zero);
}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelFlash::computeDt() {
  CH_assert(allDefined());
  // DEV CD: We compute dt in FLASH's Driver_computeDt.
  // computeDt method should not be called.
  MayDay::Error("AMRLevelFlash::computeDt: Should not be called by FLASH");
  return 0.0;
}

// Compute dt using initial data
Real AMRLevelFlash::computeInitialDt() {
  CH_assert(allDefined());
  // DEV CD: We compute dt in FLASH's Driver_computeDt.
  // computeInitialDt method should not be called.
  MayDay::Error("AMRLevelFlash::computeInitialDt: Should not be called by FLASH");
  return 0.0;
}

const LevelData<FArrayBox>& AMRLevelFlash::getStateNew() const {
  CH_assert(allDefined());
  return m_UNew;
}

const LevelData<FArrayBox>& AMRLevelFlash::getStateOld() const {
  CH_assert(allDefined());
  return m_UOld;
}

bool AMRLevelFlash::allDefined() const {
  return isDefined() && m_paramsDefined;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelFlash::loadBalance(const Vector<Box>& a_grids) {
   // assert elimintated for restart, very hackish
  //CH_assert(allDefined());

  // Load balance and create boxlayout
  Vector<int> procMap;
  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4) {
    pout() << "AMRLevelFlash::loadBalance: procesor map: " << endl;
    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
      pout() << igrid << ": " << procMap[igrid] << "  " << endl;
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();
  return dbl;
}

// Setup menagerie of data structures
void AMRLevelFlash::levelSetup() {
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    pout() << "AMRLevelFlash::levelSetup " << m_level << endl;
  
  AMRLevelFlash* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelFlash* amrGodFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);

  if (m_hasCoarser) {
    int nRefCrse = m_coarser_level_ptr->refRatio();
    m_coarseAverage.define(m_grids, m_numStates, nRefCrse);
    m_fineInterp.define(m_grids, m_numStates, nRefCrse, m_problem_domain);

    const DisjointBoxLayout& coarserLevelDomain = amrGodCoarserPtr->m_grids;
    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    if (m_useFluxCorrection) {
      amrGodCoarserPtr->m_fluxRegister.define(
        m_grids, amrGodCoarserPtr->m_grids, m_problem_domain, amrGodCoarserPtr->m_ref_ratio, m_numStates, m_scaleFineFluxes
      );
      amrGodCoarserPtr->m_fluxRegister.setToZero();
    }
  }

  // DEV CD: Create a scratch_ctr data structure.
  CH_assert (m_meshInfo.find(SCRATCH_CTR) != m_meshInfo.end());
  if (m_meshInfo[SCRATCH_CTR].size() > 0) {
    IntVect ivGhost = m_numGhost*IntVect::Unit;
    m_scratchCtrFab.define(m_grids,m_meshInfo[SCRATCH_CTR].size(),ivGhost);
  }

  // DEV CD: Create a box lookup table so that FLASH can provide
  // an integer block ID and receive a Chombo box dataIndex.
  makeBoxLookupTable();
}

// Get the next coarser level
AMRLevelFlash* AMRLevelFlash::getCoarserLevel() const {
  CH_assert(allDefined());
  AMRLevelFlash* amrGodCoarserPtr = NULL;
  if (m_coarser_level_ptr != NULL) {
    amrGodCoarserPtr = dynamic_cast<AMRLevelFlash*>(m_coarser_level_ptr);
    if (amrGodCoarserPtr == NULL)
      MayDay::Error("AMRLevelFlash::getCoarserLevel: dynamic cast failed");
  }
  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelFlash* AMRLevelFlash::getFinerLevel() const {
  CH_assert(allDefined());
  AMRLevelFlash* amrGodFinerPtr = NULL;
  if (m_finer_level_ptr != NULL) {
    amrGodFinerPtr = dynamic_cast<AMRLevelFlash*>(m_finer_level_ptr);
    if (amrGodFinerPtr == NULL)
      MayDay::Error("AMRLevelFlash::getFinerLevel: dynamic cast failed");
  }
  return amrGodFinerPtr;
}


// Methods for FLASH that are not implemented by AMRLevelPolytropicGas
// -----------------------------------------------------------------------

// Things to do after regrid
void AMRLevelFlash::postRegrid(int a_base_level) {
  CH_assert(allDefined());
  if (s_verbosity >= 3)
    pout() << "AMRLevelFlash::postRegrid " << m_level << endl;
  // DEV CD - identical to postInitialize because I expect we want same behavior
  postInitialize();
}


// New methods introduced for FLASH.
// -----------------------------------------------------------------------

void AMRLevelFlash::zeroFluxData() {
  if (m_hasFiner && m_useFluxCorrection) {
    // Clear flux registers with next finer level
    LevelFluxRegister *finerFluxRegister = &m_fluxRegister;
    finerFluxRegister->setToZero();
  }
}


void AMRLevelFlash::putFluxData(const void *pFluxes, const Real dt, const int levelBoxID, const int axis) {
  const int idir = axis;
  const int gridStruct = CENTER;
  bool debugFlux;

#ifdef DEBUG_FLUXES
  debugFlux = true;
#else
  debugFlux = false;
#endif

  CH_assert (idir >= 0 && idir < SpaceDim);

  if (m_useFluxCorrection) {
    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    const DataIndex datInd = getDataIndex(levelBoxID, gridStruct);
    const DisjointBoxLayout& dbl = getDisjointBoxLayout(gridStruct);
    const Box& b = dbl[datInd];

    // The following is copied from FluxBox::define
    // FLASH calls putFluxData one direction at a time meaning we have
    // no need for the FluxBox abstraction.
    const Box edgeBox(surroundingNodes(b,idir));
    FArrayBox F(edgeBox, m_numStates);
    // F.setVal(0.0); DEV CD: Perhaps useful for debugging purposes later.

    // DEV CD.  These debugging prints are crazy in anything more than 1D.
    // I'm running very small simulations and checking FLASH Grid_putFluxData
    // data matches the data here.
    if (debugFlux) {
      Real *printFluxes = (Real*) pFluxes;
      const int numFluxPts = (int) edgeBox.numPts();
      for (int v=0; v<m_numStates; ++v) {
        const int s = v*numFluxPts; // start cell for variable v.
        for (int c=0; c<numFluxPts; ++c) {
          const Real flx = printFluxes[s+c];
          pout() << setw(12)
             << setprecision(3)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific)
             << flx;
        }
        pout() << endl;
      }
      pout() << endl;
    }

    // Setup an interval corresponding to the conserved variables
    const Interval UInterval(0,m_numStates-1);

    // Cast to void* to remove the const from the pFluxes pointer.
    F.linearIn((void*)pFluxes, edgeBox, UInterval);

    // Increment coarse flux register between this level and the next
    // finer level - this level is the next coarser level with respect
    // to the next finer level
    if (m_hasFiner) {
      // Recall that my flux register goes between my level and the next finer level
      LevelFluxRegister* finerFluxRegister = &m_fluxRegister;
      finerFluxRegister->incrementCoarse(F,dt,datInd,
					 UInterval,
					 UInterval,idir);
      if (s_verbosity >= 4) {
        pout() << "Incremented coarse fluxes for box " << F.box() <<
          " at level " << m_level << " (see below)" << endl;
        finerFluxRegister->poutCoarseRegisters();
      }
    }

    // Increment fine flux registers between this level and the next
    // coarser level - this level is the next finer level with respect
    // to the next coarser level
    if (m_hasCoarser) {
      AMRLevelFlash* coarserPtr = getCoarserLevel();
      // Recall that my flux register goes between my level and the next finer level
      LevelFluxRegister* coarserFluxRegister = &coarserPtr->m_fluxRegister;
      coarserFluxRegister->incrementFine(F,dt,datInd, UInterval, UInterval,idir);
      if (s_verbosity >= 4) {
        pout() << "Incremented fine fluxes for box " << F.box() <<
          " at level " << m_level << " (see below)" << endl;
        coarserFluxRegister->poutFineRegisters();
      }
    }
  }
}


void AMRLevelFlash::reflux() {
    if(m_hasFiner && m_useFluxCorrection) {
        // Reflux
        Real scale = -1.0/m_dx;

        if(!m_scaleFineFluxes)
            // This code is executed when hy_fluxRepresentation="fluxes"
            // Note: This assumes dx == dy == dz
            scale = -pow(m_dx, -SpaceDim);

        m_fluxRegister.reflux(m_UNew, scale);

        // DEV CD: Not sure whether this restriction should be here.
        // Average from finer level data
        AMRLevelFlash* amrGodFinerPtr = getFinerLevel();
        amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew, amrGodFinerPtr->m_UNew);
    }
}


void AMRLevelFlash::fillGuardCells()
{
  // If there is a coarser level interpolate undefined ghost cells
  // This affects the blocks that exist on a fine-coarse interface
  if (m_hasCoarser) {
    AMRLevelFlash* amrLevelCoarserPtr = getCoarserLevel();
    const DisjointBoxLayout& cDbl = amrLevelCoarserPtr->getDisjointBoxLayout(CENTER);
    const LevelData<FArrayBox>& cData = amrLevelCoarserPtr->getStateNew();
    const ProblemDomain& cProbDom = amrLevelCoarserPtr->problemDomain();
    const Real time_interp_coeff = 0.0; // DEV CD: Need to check 0.0.

    /* Conversation with Dan Graves:
       We have two options for filling guardcells when there is a
       coarse-fine boundary: either we use piecewise linear interpolation
       (PiecewiseLinearFillPatch) or quadratic interpolation (QuadCFInterp).
       Ideally we would use QuadCFInterp, but this only fills face
       guard cell regions and does not fill corner guard cell regions.

       Fine box guard cell regions: X is filled and 0 is not filled.
       ________________
       | 0 | X X X X X
       ----------------
       | X | 
       | X |  Fine Box
       | X |
       | X |

       This is OK for elliptic problems because we do not use corner
       guard cell data, but it is an issue for hyperbolic problems
       because they use corner guard cell data.

       The general strategy implemented here is:
       1). Use piecewise linear interpolation to fill all guard cells
       on the coarse-fine boundary including the corner guard cell data.
       This ensures corner guard cells contain filled data.
       2). Use quadratic interpolation to improve the quality of the
       data at the face guard cells.  The corner guard cell data
       is untouched.
       
       DEV: Chris thought that PARAMESH always performed quadratic
       interpolation even for corner guard cell regions, but this
       is obviously not the case. */

    PiecewiseLinearFillPatch PLFP(m_grids, cDbl, m_numStates, cProbDom, m_ref_ratio, m_numGhost);
    PLFP.fillInterp(m_UNew, cData, cData, time_interp_coeff, 0, 0, m_numStates);

    if (m_useQuadCFInterp) {
      QuadCFInterp QCFI(m_grids, &cDbl, m_dx, m_ref_ratio, m_numStates, m_problem_domain);
      QCFI.coarseFineInterp(m_UNew, cData);
    }
  }

  // Fill ghost cells from boxes at the same refinement level
  m_UNew.exchange();
}


void AMRLevelFlash::coarseAverage() {
  if (m_hasFiner) {
    // Volume weighted average from finer level data
    AMRLevelFlash* amrLevelFinerPtr = getFinerLevel();
    amrLevelFinerPtr->m_coarseAverage.averageToCoarse(m_UNew, amrLevelFinerPtr->m_UNew);
  }
}


int AMRLevelFlash::getNumBoxes(const int procID, const int gridStruct) const {
  const DisjointBoxLayout& dbl = getDisjointBoxLayout(gridStruct);  
  return dbl.numBoxes(procID);
}


const DisjointBoxLayout& AMRLevelFlash::getDisjointBoxLayout(const int gridStruct) const {
  // We have a fixed reference to the disjoint box layout, m_grids.
  // This is OK for CENTER and SCRATCH_CTR because they are created from
  // the same disjoint box layout.  This will not be OK when we add
  // additional FLASH grid data structures.
  CH_assert(allDefined());
  CH_assert(gridStruct == CENTER || gridStruct == SCRATCH_CTR);
  return m_grids;
}


DataIndex AMRLevelFlash::getDataIndex(const int levelBoxID, const int gridStruct) const {
    // CD. Translate a (levelBoxID, gridStruct) pair into a dataIndex so that
    // we can directly access the appropriate Chombo patch.
    map_data_index_t::const_iterator it = m_boxLookupTable.find(gridStruct);
    CH_assert (it != m_boxLookupTable.end());
    const Vector<DataIndex> &v = it->second;
    CH_assert (levelBoxID < v.size());

    const DataIndex datInd = v[levelBoxID];
    return datInd;
}


// The pointer to grid data is the only way that FLASH can access and alter
// any m_UNew or m_scratchCtrFab data
void * AMRLevelFlash::getDataPtr(const int levelBoxID, const int gridStruct) {
  LevelData<FArrayBox> aliasLevFab;
  switch (gridStruct) {
  case CENTER:
    CH_assert (m_UNew.isDefined());
    aliasLevelData(aliasLevFab, &m_UNew, m_UNew.interval());
    break;
  case SCRATCH_CTR:
    CH_assert (m_scratchCtrFab.isDefined());
    aliasLevelData(aliasLevFab, &m_scratchCtrFab, m_scratchCtrFab.interval());
    break;
  default:
    MayDay::Error("[GetDataPtr]: Invalid mesh type");
  }

  const DataIndex datInd = getDataIndex(levelBoxID, gridStruct);
  const FArrayBox& localBlkData = aliasLevFab[datInd];
  return (void*) localBlkData.dataPtr(0);
}


box_info_t AMRLevelFlash::getBoxInfo(const int levelBoxID, const int gridStruct) const {
  // We have a fixed reference to the disjoint box layout, m_grids.
  // This is OK for CENTER and SCRATCH_CTR because they are created from
  // the same disjoint box layout.  This will not be OK when we add
  // additional FLASH grid data structures.
  CH_assert (gridStruct == CENTER || gridStruct == SCRATCH_CTR);

  const DataIndex datInd = getDataIndex(levelBoxID, gridStruct);
  box_info_t box_info;

  // Zero all data for the case where NDIM < MDIM
  for(int i=0; i<MDIM; ++i) {
    box_info.deltas[i] = 0.0;
    box_info.bsize[i] = 0.0;
    box_info.coord[i] = 0.0;
    box_info.lowBndBox[i] = 0.0;
    box_info.highBndBox[i] = 0.0;
    box_info.lowLimits[i] = 0;
    box_info.highLimits[i] = 0;
    box_info.guardcells[i] = 0;
    box_info.corner[i] = 0; //Needs knowledge of max refinement level
    box_info.stride[i] = 0; //Needs knowledge of max refinement level
    box_info.isNextToLowDomain[i] = FLASH_TRUE;
    box_info.isNextToHighDomain[i] = FLASH_TRUE;
  }

  // Some information is specific to the level rather than the box
  box_info.lrefine = m_level;
  for(int i=0; i<SpaceDim; ++i)
    box_info.deltas[i] = m_dx;

  // Other information is specific to the box
  const DisjointBoxLayout& dbl = getDisjointBoxLayout(gridStruct);
  const Box& b = dbl[datInd];
  const IntVect& small = b.smallEnd(); // e.g. (32,0)
  const IntVect& big = b.bigEnd(); // e.g. (63,31)
  const IntVect& size = b.size(); // e.g. (32,32)

  // I want to access the ProblemDomain index information.
  // This requires a box and its smallEnd and BigEnd methods.
  const Box& domainBox = m_problem_domain.domainBox();
  const IntVect& small_domain = domainBox.smallEnd();
  const IntVect& big_domain = domainBox.bigEnd();

  // Calculate the corner ID for the block of interest.
  // This should probably be moved to chombo_adaptive_grid class.
  int maxLevel = m_flashAMRInfo.maxRefineLevel;
  int refFactor = 1;
  for (int level = m_level; level < maxLevel; ++level)
    refFactor *= m_flashAMRInfo.refRatio;
  const Box br = refine(b, refFactor);
  const IntVect& smallr = br.smallEnd();

  // Note that all indexing is kept zero-based.
  // We convert to unit-based in chombo_f_c_api.
  for(int i=0; i<SpaceDim; ++i) {
    box_info.lowLimits[i] = small[i];
    box_info.highLimits[i] = big[i];
    box_info.guardcells[i] = m_numGhost;
    box_info.corner[i] = smallr[i];
    box_info.stride[i] = refFactor;
    box_info.bsize[i] = size[i] * m_dx;
    box_info.lowBndBox[i] = small[i] * m_dx;
    box_info.highBndBox[i] = (big[i] + 1) * m_dx;
    box_info.coord[i] = box_info.lowBndBox[i] +
      ((box_info.highBndBox[i] - box_info.lowBndBox[i]) / 2.0);
    box_info.isNextToLowDomain[i] =
      (small[i] == small_domain[i] ? FLASH_TRUE : FLASH_FALSE);
    box_info.isNextToHighDomain[i] =
      (big[i] == big_domain[i] ? FLASH_TRUE : FLASH_FALSE);
  }
  return box_info;
}


void AMRLevelFlash::makeBoxLookupTable() {
    Vector<DataIndex> v;
    int meshDefineVal, meshNumVars;

    m_boxLookupTable.clear();

    for(int i=0; i<MAX_GRID_DATA_STRUCT_TMP; ++i) {
        meshDefineVal = m_flashAMRInfo.meshTypes[i];
        meshNumVars = m_flashAMRInfo.meshNumVars[i];

        if(meshNumVars > 0) {
            const DisjointBoxLayout &dbl = getDisjointBoxLayout(meshDefineVal);
            DataIterator dit = dbl.dataIterator();
            for(dit.begin(); dit.ok(); ++dit)
                v.push_back(dit());
        }

        m_boxLookupTable[meshDefineVal] = v;
        v.clear();
    }
}
