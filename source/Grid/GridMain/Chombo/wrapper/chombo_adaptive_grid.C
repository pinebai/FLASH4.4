#include "chombo_adaptive_grid.h"

// Rule 1: All indices are zero-based.

// Rule 2: All vectors passed from or returned to FLASH have type std::vector.
// This will allow us to use the same interfaces for e.g. Samrai.

// Rule 3: All vectors internal to our Chombo object have type
// Vector so that Chombo can track memory.

// Note: We can find out MPI rank and communicator size using
// numProc() and procID() functions.

Chombo_Adaptive_Grid::Chombo_Adaptive_Grid()
{
  m_top_level = -1;
  m_max_level = -1;
  m_finest_level = -1;
  m_max_grid_size = 0;
  m_max_base_grid_size = m_max_grid_size;
  m_blockFactor = 4;
  m_amrlevels.resize(0);
  m_verbosity = 4;
  AMRLevel::verbosity(m_verbosity);
}


Chombo_Adaptive_Grid::~Chombo_Adaptive_Grid()
{
  for (int lev = 0; lev < m_amrlevels.size(); ++lev)
  {
    if (m_amrlevels[lev] != NULL)
    {
      delete m_amrlevels[lev];
      m_amrlevels[lev] = NULL;
    }
  }
}

void Chombo_Adaptive_Grid::Define(
    const flash_amr_info_t& flashAMRInfo,
	const mesh_info_t& meshInfo,
    bool restart
) {
    const IntVect domainVect(D_DECL6(
        flashAMRInfo.baseDomainSize[0],
        flashAMRInfo.baseDomainSize[1],
        flashAMRInfo.baseDomainSize[2],
        0,0,0));
    const IntVect ghostVect(D_DECL6(
        flashAMRInfo.guardCells[0],
        flashAMRInfo.guardCells[1],
        flashAMRInfo.guardCells[2],
        0,0,0));
    ProblemDomain computationalDomain;
    const int refRatio = flashAMRInfo.refRatio;

    m_max_level = flashAMRInfo.maxRefineLevel;
    CH_assert(m_max_level >= 0);

    if(!restart) m_finest_level = 0; // set to zero to start
    m_top_level = 0; // set to zero to start
    m_blockFactor = flashAMRInfo.BRMeshRefineBlockFactor;
    m_max_grid_size = flashAMRInfo.maxBlockSize[0];
    m_max_base_grid_size = m_max_grid_size;
    m_verbosity = flashAMRInfo.verbosity;
    m_restrictBeforeGhostExchange = FLASH_TRUE == flashAMRInfo.restrictBeforeGhostExchange;

    m_ref_ratios.resize(0);
    m_ref_ratios.resize(m_max_level+1, refRatio);

    // Assert that the requested mesh is sensible 
    for(int d=0; d<SpaceDim; ++d) {
        CH_assert(flashAMRInfo.highDomain[d] - flashAMRInfo.lowDomain[d] >= 0.0);
        CH_assert(flashAMRInfo.maxBlockSize[d] >= 0);
        CH_assert(flashAMRInfo.maxBlockSize[0] == flashAMRInfo.maxBlockSize[d]);
        CH_assert(flashAMRInfo.baseDomainSize[d] > 0);
        CH_assert(flashAMRInfo.guardCells[d] > 0);
        CH_assert(flashAMRInfo.guardCells[0] == flashAMRInfo.guardCells[d]);
    }
    
    // Zero based coordinates so subtract unit vector from the total domain size
    computationalDomain.define(Box(IntVect::Zero, (domainVect - IntVect::Unit)));
    for(int d=0; d<SpaceDim; ++d)
        // Specify whether each dimension is periodic
        computationalDomain.setPeriodic(d, flashAMRInfo.domainBC[d] == PERIODIC);

    // Setup the FLASH factory object which produces AMR levels.
    //m_amrFlashFactory.define(flashAMRInfo, meshInfo);

    { // create the level stack.  if this was a restart, then the levels have already been created, they just need defining
        while(m_max_level+1 < m_amrlevels.size()) {
            delete m_amrlevels.back();
            m_amrlevels.resize(m_amrlevels.size()-1);
        }
        m_amrlevels.resize(m_max_level+1, NULL);
        
        ProblemDomain ref_domain = computationalDomain;
        for(int l=0; l < m_amrlevels.size(); l++) {
            AMRLevelFlash *lev;
            if(NULL == (lev = dynamic_cast<AMRLevelFlash*>(m_amrlevels[l])))
                m_amrlevels[l] = lev = new AMRLevelFlash(l);
            lev->defineParams(flashAMRInfo, meshInfo, restart);
            lev->define(l == 0 ? NULL : m_amrlevels[l-1], ref_domain, l, m_ref_ratios[l]);
            if(l > 0) m_amrlevels[l-1]->finerLevelPtr(lev);
            ref_domain.refine(m_ref_ratios[l]);
        }

        for(int l=0; l < m_amrlevels.size(); l++) {
            if(restart) {
	            dynamic_cast<AMRLevelFlash*>(m_amrlevels[l])->restartGrid();
            }
        }
    }
    
    // build BRMeshRefine object.
    Real defaultFillRatio = flashAMRInfo.BRMeshRefineFillRatio;
    int defaultBufferSize = flashAMRInfo.BRMeshRefineBufferSize;
    m_mesh_refine.define(computationalDomain, m_ref_ratios, defaultFillRatio, m_blockFactor, defaultBufferSize, m_max_grid_size);

    // All processors construct the same metadata about the global domain
    if(!restart) {
        m_old_grids.resize(1);
        makeBaseLevelMesh(m_old_grids[0]);
    }
    else {
        m_old_grids.resize(m_amrlevels.size());
        for(int l=0; l < m_amrlevels.size(); l++)
            m_old_grids[l] = m_amrlevels[l]->boxes();
        UpdateBoxLevelInfo();
    }
    // m_max_level can be zero: m_old_tags is not referenced in this case.
    m_old_tags.resize(m_max_level);
}

void Chombo_Adaptive_Grid::ReadCheckpoint(
    const char *filename,
    double &simTime,
    double &dt,
    HDF5HeaderData &scalars,
    HDF5HeaderData &runparms) {
    
    HDF5Handle h(filename, HDF5Handle::OPEN_RDONLY);
    
    h.pushGroup("Flash_scalars");
    scalars.readFromFile(h);
    h.popGroup();
    
    h.pushGroup("Flash_runparms");
    runparms.readFromFile(h);
    h.popGroup();
    
    {
        HDF5HeaderData hdr;
        hdr.readFromFile(h);
        if(!hdr.m_string.count("filetype") || hdr.m_string["filetype"] != "VanillaAMRFileType")
            MayDay::Error("Checkpoint file invalid.");
        m_finest_level = hdr.m_int["num_levels"] - 1;
    }
    
    for(int l=0; l <= m_finest_level; l++) {
        AMRLevelFlash *lev = new AMRLevelFlash(l);
        lev->readCheckpointHeader(h);
        lev->readCheckpointLevel(h);
        m_amrlevels.push_back(lev);
    }
    
    h.close();
}
void Chombo_Adaptive_Grid::WriteCheckpoint(
    const char filename[],
    double simTime,
    double dt,
    const HDF5HeaderData &scalars,
    const HDF5HeaderData &runparms) const {
#ifdef CH_USE_HDF5
  HDF5Handle handle(filename, HDF5Handle::CREATE);

  // Add a header containing extra metadata.  The "filetype"
  // HDF5 attribute does not seem to be required, but the
  // "num_levels" HDF5 attribute is definitely required - visit 
  // uses "num_levels" to recognize Chombo file format.
  HDF5HeaderData header;
  string filedescriptor("VanillaAMRFileType");
  header.m_string["filetype"] = filedescriptor;
  header.m_int["num_levels"] = m_finest_level + 1;
  header.writeToFile(handle);
  
  handle.pushGroup("Flash_scalars");
  scalars.writeToFile(handle);
  handle.popGroup();
  
  handle.pushGroup("Flash_runparms");
  runparms.writeToFile(handle);
  handle.popGroup();

  for (int ilev = 0; ilev <= m_finest_level; ilev++) {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->time(simTime); // FLASH sets time
    amrLevelFlashPtr->dt(dt); // FLASH sets dt
    amrLevelFlashPtr->writeCheckpointHeader(handle);
    amrLevelFlashPtr->writeCheckpointLevel(handle);
  }

  handle.close();
#endif
}


// Method depends on class level data:
// m_max_base_grid_size, m_blockFactor, m_amrlevels
void Chombo_Adaptive_Grid::makeBaseLevelMesh(Vector<Box>& a_grids) const
{
  if (m_max_base_grid_size == 0)
  {
    // define base level to be single grid
    a_grids.resize(1);
    a_grids[0] = m_amrlevels[0]->problemDomain().domainBox();
  }
  else
  {
    if (m_max_base_grid_size < m_blockFactor)
    {
      MayDay::Abort("Base grid size must be greater than blocking factor");
    }

    // chop base level up into grids of no more than m_max_grid_size on a side
    // first coarsen to enforce the blocking factor
    ProblemDomain problem_domain = coarsen(m_amrlevels[0]->problemDomain(),
					   m_blockFactor);
    if (refine(problem_domain, m_blockFactor) != m_amrlevels[0]->problemDomain())
    {
      MayDay::Error("level 0 problem domain not coarsenable by blocking factor");
    }
    int max_grid_size = m_max_base_grid_size / m_blockFactor;
    Tuple<Vector<int>,SpaceDim> box_sizes;
    IntVect num_grids;
    IntVect base_size;

    for (int d = 0; d < SpaceDim; ++d)
    {
      int num_div = 1;
      int domain_size = problem_domain.domainBox().size(d);
      while (num_div * max_grid_size < domain_size)
      {
	++num_div;
      }

      // int(x/y) +(x%y)?1:0 is integer division with rounding upwards, for x,y>0
      base_size[d] = int(domain_size/num_div) +((domain_size%num_div) ? 1 : 0);
      box_sizes[d].resize(num_div, base_size[d]);
      box_sizes[d][num_div-1] = domain_size -(num_div - 1) * base_size[d];
      num_grids[d] = num_div;
    }

    Box b(IntVect::Zero,num_grids - IntVect::Unit);
    const IntVect& domain_hi = problem_domain.domainBox().bigEnd();
    const IntVect& domain_lo = problem_domain.domainBox().smallEnd();

    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      IntVect lo = domain_lo + iv * base_size;
      IntVect hi = min(lo + base_size - IntVect::Unit, domain_hi);
      Box grid(lo,hi);

      a_grids.push_back(refine(grid, m_blockFactor) );
    }
  }
}


// If False we need to call BRMeshRefine::regrid again.
bool Chombo_Adaptive_Grid::IsInitialRefinementDone()
  const
{
  bool initialRefinementDone;
  initialRefinementDone = !(m_top_level < Min(m_finest_level+1, m_max_level));
  if (m_verbosity >= 2)
    pout() << "Chombo_Adaptive_Grid::IsInitialRefinementDone() " <<
      initialRefinementDone << endl;
  return initialRefinementDone;
}


// FLASH executes the do while loop for level generation.
void Chombo_Adaptive_Grid::BuildInitialGrid()
{
  if (m_verbosity >= 2)
    pout() << "Chombo_Adaptive_Grid::BuildInitialGrid()" << endl;

  if (IsInitialRefinementDone())
  {
    // The mesh is now appropriately refined.  Call AMRLevel::initialGrid
    // one last time for each level to load balance the blocks and
    // allocate the solution data arrays.
    if (m_finest_level == 0)
    {
      m_new_grids.resize(1);
      m_new_grids[0] = m_old_grids[0];
    }

    for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialGrid(m_new_grids[level]);
    }

    if (m_verbosity >= 2)
      PrintBoxes("BuildInitialGrid (final boxes)", m_new_grids);
  }
  else
  {
    for (int level = 0; level <= m_top_level; ++level)
    {
      m_amrlevels[level]->initialGrid(m_old_grids[level]);
    }

    if (m_verbosity >= 2)
      PrintBoxes("BuildInitialGrid (intermediate boxes)", m_old_grids);
  }

  // Update FLASH block to Chombo box mapping.
  UpdateBoxLevelInfo();
}


void Chombo_Adaptive_Grid::RefineInitialGrid()
{
  if (m_verbosity >= 2)
    pout() << "Chombo_Adaptive_Grid::RefineInitialGrid()" << endl;

  // The FLASH subroutine Grid_markRefineDerefine has marked the
  // tagc mesh variable of each cell in all blocks with a 0.0 or 1.0.
  // A 1.0 means the cell is marked for refinement.  I use the
  // tagCellsInit function to construct an IntVectSet by iterating
  // over each cell and checking whether it is marked with a 1.0.
  m_tags.clear();
  m_tags.resize(m_finest_level+1);
  
  for (int lev=0; lev<=m_finest_level; ++lev)
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[lev]);
    amrLevelFlashPtr->tagCellsInit(m_tags[lev]);
  }

  // m_new_grids is resized and filled in regrid function
  // m_tags is also modified.  All other args are constant
  m_finest_level = m_mesh_refine.regrid(m_new_grids, m_tags, 0, m_top_level, m_old_grids);

  // do this only if a new level was generated
  if (m_finest_level > m_top_level)
  {
    m_old_grids = m_new_grids;
    // copy current tags to old_tags
    for (int lev = 0; lev <= m_top_level; lev++)
    {
      m_old_tags[lev] = m_tags[lev];
    }

    if (m_verbosity >= 2)
      PrintBoxes("RefineInitialGrid (new level added)", m_old_grids);
  }

  // Note that m_new_grids is a Vector<Vector<Box> > so does not have a
  // process mapping yet.  Do not call UpdateBoxLevelInfo() yet.

  // The variable m_top_level is used in an identical way to the loop index
  // named top_level in AMR::initialGrid.  The name is slightly
  // misleading because it is not actually the current top level, but I will
  // leave the name alone so my code looks more similar to AMR::initialGrid.
  ++m_top_level;
}


void Chombo_Adaptive_Grid::FinalizeInitialGrid()
{
  if (m_verbosity >= 2)
    pout() << "Chombo_Adaptive_Grid::FinalizeInitialGrid()" << endl;

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
  {
    m_amrlevels[level]->initialGrid(Vector<Box>());
  }

  // call post-initialize once all the levels have been defined
  for (int level = m_finest_level; level >= 0; --level)
  {
    m_amrlevels[level]->postInitialize();
  }
}


// generate new grid hierarchy
void Chombo_Adaptive_Grid::Regrid(const int a_base_level)
{
  CH_TIME("Chombo_Adaptive_Grid::Regrid");

  if (m_verbosity >= 2)
  {
    pout() << "Chombo_Adaptive_Grid::Regrid(" << a_base_level << ")" << endl;
  }

  int top_level = Min(m_finest_level, m_max_level - 1);

  // DEV CD - Most of this code is copied from AMR::regrid.
  // I added the line below to ensure that top_level is >= 0.
  // This stops problems when m_finest_level = m_max_level = 0.
  // I also added assertions about m_max_level and m_finest_level.
  top_level = Max(0, top_level);
  CH_assert(m_max_level >= 0);
  CH_assert(m_finest_level >= 0);
  if(m_verbosity >= 2)
    pout() << "Chombo_Adaptive_Grid::Regrid:"
        << "  top_level = " << top_level << endl
        << "  m_max_level = " << m_max_level << endl
        << "  m_finest_level = " << m_finest_level << endl;

  Vector<Vector<Box> > old_grids(top_level+1);
  Vector<Vector<Box> > new_grids;
  Vector<IntVectSet> tags(top_level+1);
  Vector<ProblemDomain> problem_domains(top_level+1);

  for (int level = a_base_level; level <= top_level; ++level) {
    m_amrlevels[level]->tagCells(tags[level]);
    old_grids[level] = m_amrlevels[level]->boxes();
    problem_domains[level] = m_amrlevels[level]->problemDomain();

    if(m_verbosity >= 2)
        pout() << "Chombo_Adaptive_Grid::Regrid: problem domain[" << level << "]: " << problem_domains[level] << endl;

    if(m_verbosity >= 4) {
      pout() << "Chombo_Adaptive_Grid::Regrid: old_grids[" << level << "]: " << endl;
      for(int i = 0; i< old_grids[level].size(); ++i)
        pout() << "  " << i << ": " << old_grids[level][i] << endl;
      if(m_verbosity >= 4)
        pout() << "Chombo_Adaptive_Grid::Regrid: tags[" << level << "]: " << tags[level] << endl;
    }
  }

  int new_finest_level = m_mesh_refine.regrid(new_grids, tags, a_base_level, top_level, old_grids);
  //can only add one level at a time
  new_finest_level = Min(m_finest_level+1, new_finest_level);

  if((m_finest_level != new_finest_level) && (m_verbosity >= 2))
    pout() << "finest level changes here from " << m_finest_level << " to " << new_finest_level << endl;

  //allow for levels to change
  m_finest_level = Min(new_finest_level, m_max_level);

    if(m_verbosity >= 4) {
        if(new_grids.size() == 0)
            pout() << "No new_grids" << endl << endl;
        else {
            for(int level = a_base_level; level <= m_finest_level; ++level) {
                pout() << "new_grids[" << level << "]: " << endl;
                for(int i = 0; i< new_grids[level].size(); ++i)
                    pout() << "  " << i << ": " << new_grids[level][i] << endl;
                pout() << endl;
            }
        }
    }

  // before regridding (but after tagging) allow for pre-regridding ops
  // (tjl 8/23/06 - this needed for mapped grid computations and a default
  // implementation is provided in AMRLevel to preserve existing codes)
  for (int level = m_finest_level; level >= a_base_level; --level)
    m_amrlevels[level]->preRegrid(a_base_level, new_grids);

  for (int level = a_base_level + 1; level <= m_finest_level; ++level)
    m_amrlevels[level]->regrid(new_grids[level]);

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    m_amrlevels[level]->regrid(Vector<Box>());

  // now that the new hierarchy is defined, do post-regridding ops
  // (dfm 8/26/05 -- call postRegrid on base_level as well, to
  // cover the case where all levels finer than base_level are removed)
  for (int level = m_finest_level; level >= a_base_level; --level)
    m_amrlevels[level]->postRegrid(a_base_level);

  // Update FLASH block to Chombo box mapping.
  UpdateBoxLevelInfo();

  if (m_verbosity >= 2) {
    PrintBoxes("Regrid event (new grid shown)", new_grids);
    const size_t nBlk = m_boxLevelInfo.size(); // total boxes on myPE
    for (size_t i=0; i<nBlk; ++i)
      pout() << "Local box " << i <<
	" has level " <<  m_boxLevelInfo[i].level <<
	" and level box ID " << m_boxLevelInfo[i].levelBoxID << endl;
  }

  // DEV CD: Eventually we want to know when the mesh actually changes.
  // This does not seem to be trivial.  Some possibilities:
  //
  // 1. Compare the boxes in new_grids and old_grids.  We need to be
  // careful because old_grids is load balanced and new_grids is not.
  //
  // 2. Add a preRegrid method in AMRLevelFLash which stores a deep copy
  // of m_grids (a DisjointBoxLayout) before it is modified in regrid.
  // Both these DisjointBoxLayouts will be load balanced and we can extract
  // the boxes using the boxArray method of BoxLayout.  One concern is that
  // there is a comment above deepCopy method in BoxLayout.H that says:
  // "There is no assurance that the order in which this BoxLayout
  // is indexed corresponds to the indexing of <i>a_source</i>."
  //
  // I think for 1 and 2 we have to sort the boxes before we can check
  // whether the boxes in the mesh are the same.
}


std::vector<int> Chombo_Adaptive_Grid::GetBlockIDs(const int blkType,
						   const int level)
  const
{
  // The indices of vBlk are used to access a cached box lookup table.
  // It gives FLASH fast access to Chombo boxes.  The indices give
  // access to all the boxes that match blkType (and level) criteria.
  const int gridStruct = CENTER;
  std::vector<int> vBlk;
  const size_t nBlk = m_boxLevelInfo.size(); // total boxes on myPE

  if (blkType == ALL_BLKS || blkType == LEAF)
  {
    // Return all indices - the most common situation
    for (size_t i=0; i<nBlk; ++i)
    {
      vBlk.push_back(i);
    }
  }
  else if (blkType == REFINEMENT)
  {
    // Return all indices for boxes at a chosen level
    CH_assert (level >= 0 && level <= m_max_level);
    for (size_t i=0; i<nBlk; ++i)
    {
      if (level == m_boxLevelInfo[i].level)
      {
	vBlk.push_back(i);
      }
    }
  }
  else
  {
    // Return all indices for boxes next to chosen boundaries
    box_info_t boxInfo;
    int s, e;

    switch (blkType) {
    case IBDRY_BLKS:
      s = e = 0; // tests only IAXIS
      break;
    case JBDRY_BLKS:
      s = e = 1; // tests only JAXIS
      break;
    case KBDRY_BLKS:
      s = e = 2; // tests only KAXIS
      break;
    case ANY_BDRY_BLKS:
      s = 0; e = 2; // tests IAXIS, JAXIS and KAXIS
      break;
    default:
      MayDay::Error("[GetBlockIDs]: Invalid boundary type");
    }

    for (size_t i=0; i<nBlk; ++i)
    {
      boxInfo = GetBoxInfo(i, gridStruct);

      // Test whether block i is touching dimension j boundary
      for (int j=s; j<=e; ++j)
      {
	if (boxInfo.isNextToLowDomain[j] == FLASH_TRUE ||
	    boxInfo.isNextToHighDomain[j] == FLASH_TRUE)
	{
	  vBlk.push_back(i);
	  break; // So we do not add the same block multiple times
	}
      }
    }
  }
  return vBlk;
}


std::vector<double> Chombo_Adaptive_Grid::GetCellCoords(const int blkID,
							const int axis,
							const int edge,
							const int guardcell)
  const
{
  const int gridStruct = CENTER;
  std::vector<double> vCoords;
  box_info_t box_info;
  double lowBnd, delta, halfDelta, coords;
  int cell, numCells, numGuardCells;

  CH_assert(axis >= 0 && axis < SpaceDim);
  CH_assert(edge == LEFT_EDGE || edge == CENTER || 
	    edge == RIGHT_EDGE || edge == FACES);
  CH_assert(guardcell == 0 || guardcell == 1);

  // Pass the hard-coded gridStruct = CENTER.  Use the argument named
  // edge to find appropriate coordinates for CENTER and FACES.
  box_info = GetBoxInfo(blkID, gridStruct);
  cell = 0;
  numCells = box_info.highLimits[axis] - box_info.lowLimits[axis] + 1;
  lowBnd = box_info.lowBndBox[axis];
  delta = box_info.deltas[axis];
  halfDelta = delta / 2.0;
  numGuardCells = box_info.guardcells[axis];

  if (guardcell == 1)
  {
    cell -= numGuardCells;
    numCells += (2 * numGuardCells);
  }

  for (int i=0; i<numCells; ++i)
  {
    if (edge == LEFT_EDGE || edge == FACES)
    {
      coords = lowBnd + (cell * delta);
    }
    if (edge == CENTER)
    {
      coords = lowBnd + (cell * delta) + halfDelta;
    }
    ++cell;
    if (edge == RIGHT_EDGE)
    {
      coords = lowBnd + (cell * delta);
    }
    vCoords.push_back(coords);
  }
  if (edge == FACES)
  {
    coords = lowBnd + (cell * delta);
    vCoords.push_back(coords);
  }

  return vCoords;
}


box_info_t Chombo_Adaptive_Grid::GetBoxInfo(const int blkID, const int gridStruct) const {
  box_level_info_t boxLevelInfo = m_boxLevelInfo[blkID];
  AMRLevelFlash* amrLevelFlashPtr =
    dynamic_cast<AMRLevelFlash*>(m_amrlevels[boxLevelInfo.level]);
  box_info_t boxInfo =
    amrLevelFlashPtr->getBoxInfo(boxLevelInfo.levelBoxID,gridStruct);
  return boxInfo;
}


// This function should be called whenever the configuration of the mesh
// changes.  This means within Chombo_Adaptive_Grid::BuildInitialGrid
// and Chombo_Adaptive_Grid::Regrid.

// UpdateBoxLevelInfo can only be called after AMRLevel::initialGrid
// because AMRLevel::initialGrid creates the load balanced DisjointBoxLayout
// that UpdateBoxLevelInfo queries.  Before AMRLevel::initialGrid we only
// have a vector of boxes with no processor mapping.

void Chombo_Adaptive_Grid::UpdateBoxLevelInfo()
{
  const int gridStruct = CENTER;
  // Fixed gridStruct because the number of blocks is the same for each
  // grid data structure.
  m_boxLevelInfo.clear();

  for (int ilev=0; ilev<=m_finest_level; ++ilev)
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    const int numBoxes = amrLevelFlashPtr->getNumBoxes(procID(), gridStruct);

    for (int i=0; i<numBoxes; ++i)
    {
      box_level_info_t boxLevelInfo;
      boxLevelInfo.level = ilev;
      boxLevelInfo.levelBoxID = i;

      // This data structure will be used as follows:
      // 1). FLASH requests box ID 24.
      // 2). (level, levelBoxID) = m_boxLevelInfo[24].
      // 3). amrLevelFlash[level].getBlkPtr[levelBoxID].
      m_boxLevelInfo.push_back(boxLevelInfo);
    }
  }
}


void * Chombo_Adaptive_Grid::GetDataPtr(const int blkID,
					const int gridStruct)
{
  box_level_info_t boxLevelInfo = m_boxLevelInfo[blkID];
  AMRLevelFlash* amrLevelFlashPtr = 
    dynamic_cast<AMRLevelFlash*>(m_amrlevels[boxLevelInfo.level]);
  return amrLevelFlashPtr->getDataPtr(boxLevelInfo.levelBoxID, gridStruct);
}


// Fill ghost cells with valid data.
// DEV CD: I have looked at execPolytropic Chombo application and find that
// ghost cells are exchanged in LevelGodunov.step and
// AMRLevelPolytropicGas.tagCells:
// AMR.timestep -> AMRLevelPolytropicGas.advance -> LevelGodunov.step
// AMR.regrid AND AMR.initialGrid -> AMRLevelPolytropicGas.tagCells
// The levels are iterated from coarse to fine and so I do the same here.
void Chombo_Adaptive_Grid::FillGuardCells()
{
  if (m_restrictBeforeGhostExchange)
  {
    // DEV CD: The m_restrictBeforeGhostExchange option is a hack that
    // allows us to make guard cell fills more similar to Paramesh.  This
    // should not be left in the code.  I start from m_finest_level-1
    // so that the same code works during initial grid generation.
    // m_hasFiner is always True during initial grid generation because
    // all levels are allocated first.  Subtracting 1 does not miss the
    // finest level because of the way code is structured in coarseAverage.
    for (int ilev=m_finest_level-1; ilev>=0; ilev--)  //down
    {
      AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
      amrLevelFlashPtr->coarseAverage();
    }
  }

  for(int ilev=0; ilev<=m_finest_level; ilev++)
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->fillGuardCells();
  }
}


// Average data from fine levels to coarse levels after
// each FLASH data update!!!  This means FLASH client code needs to 
// be modified to do this averaging after each update.
// DEV CD: A coarse average happens in AMRLevelPolytropicGas.postTimeStep 
// and AMRLevelPolytropicGas.postInitialize:
// AMR.initialGrid -> AMRLevelPolytropicGas.postInitialize
// AMR.timeStep -> AMRLevelPolytropicGas.postTimeStep
// The levels are iterated from fine to coarse and so I do the same here.
void Chombo_Adaptive_Grid::AverageLevelData()
{
  for (int ilev=m_finest_level; ilev>=0; ilev--)  //down
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->coarseAverage();
  }
}


void Chombo_Adaptive_Grid::ZeroFluxData()
{
  for(int ilev=0; ilev<=m_finest_level; ilev++)
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->zeroFluxData();
  }
}


void Chombo_Adaptive_Grid::PutFluxData(const void *pFluxes,
				       const Real dt,
				       const int blkID,
				       const int axis)
{
  box_level_info_t boxLevelInfo = m_boxLevelInfo[blkID];
  AMRLevelFlash* amrLevelFlashPtr = 
    dynamic_cast<AMRLevelFlash*>(m_amrlevels[boxLevelInfo.level]);
  amrLevelFlashPtr->putFluxData(pFluxes,
				dt,
				boxLevelInfo.levelBoxID,
				axis);
}


// DEV CD: The reflux in Chombo execPolytropic application happens in
// AMRLevelPolytropicGas.postTimeStep.  This is called from AMR.timeStep
// and the iteration is fine to coarse.
void Chombo_Adaptive_Grid::Reflux()
{
  for (int ilev=m_finest_level; ilev>=0; ilev--)  //down
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->reflux();
  }
}


// DEV CD: Future software engineering task is to create a higher-order
// function that takes a pointer reflux and postTimeStep method as
// an argument.  This will prevent copying and pasting the for loop and
// cast to AMRLevelFlash each time.
void Chombo_Adaptive_Grid::PostTimeStep()
{
  for (int ilev=m_finest_level; ilev>=0; ilev--)  //down
  {
    AMRLevelFlash* amrLevelFlashPtr = dynamic_cast<AMRLevelFlash*>(m_amrlevels[ilev]);
    amrLevelFlashPtr->postTimeStep();
  }
}


// Flux conservation:
// increment coarse on the coarse side.
// increment fine on the fine side.
// reflux .... Grid_conserveFlux
// flux_density in paramesh Flash.h.


// unit aspect ratio - dx (no vector of cell spaces)
// RealVect store dx.  ignore anisotropic problems

// A simple utility function that prints the hierarchy of boxes to file.
void Chombo_Adaptive_Grid::PrintBoxes(const string msg,
				      const Vector<Vector<Box> >& vvb)
  const
{
  pout() << endl << msg << endl;
  for (size_t i=0; i<vvb.size(); ++i)
  {
    const Vector<Box>& vb = vvb[i];
    for (size_t j=0; j<vb.size(); ++j)
    {
      pout() << " Level " << i
	     << " Box " << j
	     << " : " << vb[j] << endl;
    }
  }
  pout() << endl;
}
