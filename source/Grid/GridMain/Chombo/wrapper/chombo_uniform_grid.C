#include "chombo_uniform_grid.h"

// Rule 1: All indices are zero-based.

// Rule 2: All vectors passed from or returned to FLASH have type std::vector.
// This will allow us to use the same interfaces for e.g. Samrai.

// Rule 3: All vectors internal to our Chombo object have type
// Vector so that Chombo can track memory.

// Note: We can find out MPI rank and communicator size using
// numProc() and procID() functions.

Chombo_Uniform_Grid::Chombo_Uniform_Grid()
{
  for (int i=0; i<SpaceDim; ++i)
  {
    m_deltas[i] = 0.0;
    m_ghostVect[i] = 0;
  }
  m_verbosity = 0;
}


Chombo_Uniform_Grid::~Chombo_Uniform_Grid()
{
}


void Chombo_Uniform_Grid::Define(
    const flash_ug_info_t& flashUGInfo,
	const mesh_info_t& meshInfo,
    bool restart) {
        
    // Computational domain
    const IntVect domainVect(D_DECL6(flashUGInfo.baseDomainSize[0],
                   flashUGInfo.baseDomainSize[1],
                   flashUGInfo.baseDomainSize[2],0,0,0));

    // Guard cells per block
    const IntVect ghostVect(D_DECL6(flashUGInfo.guardCells[0],
                  flashUGInfo.guardCells[1],
                  flashUGInfo.guardCells[2],0,0,0));

    // Single block size
    const int iBlkSize = flashUGInfo.baseDomainSize[0] / flashUGInfo.procGrid[0];
    const int jBlkSize = flashUGInfo.baseDomainSize[1] / flashUGInfo.procGrid[1];
    const int kBlkSize = flashUGInfo.baseDomainSize[2] / flashUGInfo.procGrid[2];

    Vector<Box> boxes;
    int xl, xu, yl, yu, zl, zu;

    m_flashUGInfo = flashUGInfo;
    m_meshInfo = meshInfo;
    m_ghostVect = ghostVect;
    m_verbosity = flashUGInfo.verbosity;

    if (m_verbosity >= 2)
    PrintMeshNames();


    // Assert that the requested mesh is sensible 
    for (int d=0; d<SpaceDim; ++d) {
        CH_assert ((flashUGInfo.highDomain[d] - flashUGInfo.lowDomain[d]) >= 0.0);
        CH_assert (flashUGInfo.procGrid[d] > 0);
        CH_assert (flashUGInfo.baseDomainSize[d] > 0);
        CH_assert (flashUGInfo.guardCells[d] > 0);
        CH_assert (flashUGInfo.guardCells[0] == flashUGInfo.guardCells[d]);
        CH_assert (flashUGInfo.baseDomainSize[d] % flashUGInfo.procGrid[d] == 0);
    }
    CH_assert(flashUGInfo.procGrid[0]*flashUGInfo.procGrid[1]*flashUGInfo.procGrid[2] == numProc());

    // Zero based coordinates so subtract unit vector from the total domain size
    m_computationalDomain.define(Box(IntVect::Zero, (domainVect - IntVect::Unit)));
    for (int d=0; d<SpaceDim; ++d) {
        // Specify whether each dimension is periodic
        m_computationalDomain.setPeriodic(d, (flashUGInfo.domainBC[d] == PERIODIC));
    }

    // Calculate the cell spacing
    for (int d=0; d<SpaceDim; ++d) {
        m_deltas[d] = (flashUGInfo.highDomain[d] - flashUGInfo.lowDomain[d]) / 
            m_computationalDomain.domainBox().size(d);
    }

    if(!restart) {
        // All processors construct the same metadata about the global domain
        for (int k=0; k<flashUGInfo.procGrid[2]; ++k) {
            zl = k * kBlkSize;
            zu = zl + kBlkSize - 1;

            for (int j=0; j<flashUGInfo.procGrid[1]; ++j) {
                yl = j * jBlkSize;
                yu = yl + jBlkSize - 1;

                for (int i=0; i<flashUGInfo.procGrid[0]; ++i) {
                    xl = i * iBlkSize;
                    xu = xl + iBlkSize - 1;
                    boxes.push_back(Box(IntVect(D_DECL6(xl,yl,zl,0,0,0)), IntVect(D_DECL6(xu,yu,zu,0,0,0))));
                }
            }
        }

        // One block per processor - definition of FLASH UG
        Vector<int> procIDs;
        for (int ibox=0; ibox<boxes.size(); ++ibox) {
            CH_assert (ibox < numProc());
            procIDs.push_back(ibox);
        }
        m_ugBoxLayout.define(boxes, procIDs, m_computationalDomain);
        m_ugBoxLayout.close();
    }
    else {
        m_ugBoxLayout = m_ctrFab.disjointBoxLayout();
        m_ghostVect = m_ctrFab.ghostVect();
    }
    
    // Ensure that we have information about the following meshes:
    CH_assert (m_meshInfo.find(CENTER) != m_meshInfo.end());
    CH_assert (m_meshInfo.find(SCRATCH_CTR) != m_meshInfo.end());

    // Add guard cells to the box layout
    if (!restart && m_meshInfo[CENTER].size() > 0) {
        m_ctrFab.define(m_ugBoxLayout, m_meshInfo[CENTER].size(), m_ghostVect);
    }
    if (m_meshInfo[SCRATCH_CTR].size() > 0) {
        m_scratchCtrFab.define(m_ugBoxLayout, m_meshInfo[SCRATCH_CTR].size(), m_ghostVect);
    }

    // Since metadata about all boxes is stored on all processors, we could
    // print the global grid decomposition using just a single processor
    if (m_verbosity >= 1) {
        pout() << "\nThe UG global domain is distributed as follows:\n" << m_ugBoxLayout << endl;
    }
}


template <class T>
static int readLevelRedefine(HDF5Handle&   a_handle,
              const int&    a_level,
              LevelData<T>& a_data,
              Real& a_dx,
              Real& a_dt,
              Real& a_time,
              Box&  a_domain,
              int&  a_refRatio,
              const Interval&   a_comps,
              bool  setGhost)
{
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  //unused
  // int nComp = header.m_int["num_components"];

  int error;
  char levelName[10];
  sprintf(levelName, "level_%i",a_level);
  error = a_handle.pushGroup(levelName);
  if(error != 0) return 1;

  HDF5HeaderData meta;
  error = meta.readFromFile(a_handle);
  if(error != 0) return 2;
  a_dx       = meta.m_real["dx"];
  a_dt       = meta.m_real["dt"];
  a_time     = meta.m_real["time"];
  a_domain   = meta.m_box["prob_domain"];
  a_refRatio = meta.m_int["ref_ratio"];

  Vector<Box> boxes;
  error = read(a_handle, boxes);
  Vector<int> procIDs;
  LoadBalance(procIDs, boxes);

  DisjointBoxLayout layout(boxes, procIDs);

  layout.close();
  if(error != 0) return 3;
  error = read(a_handle, a_data, "data", layout, a_comps, true);//added true for redefine

  if(error != 0) return 4;

  a_handle.popGroup();

  return 0;
}
void Chombo_Uniform_Grid::ReadCheckpoint(
    const char *filename,
    double &simTime,
    double &dt,
    HDF5HeaderData &scalars,
    HDF5HeaderData &runparms) {
    
    int ncomps;
    
    HDF5Handle h(filename, HDF5Handle::OPEN_RDONLY);
    {
        HDF5HeaderData meta;
        meta.readFromFile(h);
        ncomps = meta.m_int["num_compoonents"];
        if(meta.m_int.find("num_levels") != meta.m_int.end() && meta.m_int["num_levels"] != 1)
            MayDay::Error("ReadChomboCheckpoint: num_levels != 1.");
    }
    
    h.pushGroup("Flash_scalars");
    scalars.readFromFile(h);
    h.popGroup();
    
    h.pushGroup("Flash_runparms");
    runparms.readFromFile(h);
    h.popGroup();
    
    int refRatio;
    Box dombox;
    readLevelRedefine(h, 0, m_ctrFab, m_deltas[0], dt, simTime, dombox, refRatio, Interval(), false);
    h.close();
}
void Chombo_Uniform_Grid::WriteCheckpoint(
    const char *filename,
    double simTime,
    double dt,
    const HDF5HeaderData &scalars,
    const HDF5HeaderData &runparms) const {
#ifdef CH_USE_HDF5
    HDF5Handle handle(filename, HDF5Handle::CREATE);

    // Write metadata
    mesh_info_t::const_iterator it = m_meshInfo.find(CENTER);
    if (it != m_meshInfo.end()) {
        HDF5HeaderData header;

        const std::vector<std::string> &v = it->second;
        int numlevels = 1;
        int nComp = v.size();
        string filedescriptor("VanillaAMRFileType");
        header.m_string ["filetype"] = filedescriptor;
        header.m_int ["num_levels"] = numlevels;
        header.m_int ["num_components"] = nComp;

        for (int ivar = 0; ivar < nComp; ivar++)
        {
          char labelChSt[80];
          sprintf(labelChSt, "component_%d", ivar);
          string label(labelChSt);
          header.m_string[label] = v[ivar];
        }

        header.writeToFile(handle);
    }
    handle.pushGroup("Flash_scalars");
    scalars.writeToFile(handle);
    handle.popGroup();

    handle.pushGroup("Flash_runparms");
    runparms.writeToFile(handle);
    handle.popGroup();

    // Write data
    Real dtLevel = dt;
    Real dxLevel = m_deltas[0];
    Real time = simTime;
    int refLevel = 1;
    int ilev = 0;
    int err = writeLevel(handle, ilev, m_ctrFab, dxLevel, dtLevel, time, m_computationalDomain.domainBox(), refLevel);
    if(err != 0) MayDay::Error("WriteChomboPlotFile: Error in writeLevel");

    handle.close();
#endif
}

std::vector<double> Chombo_Uniform_Grid::GetCellCoords(const int blkID,
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

  CH_assert(blkID == 0);
  CH_assert(axis >= 0 && axis < SpaceDim);
  CH_assert(edge == LEFT_EDGE || edge == CENTER || 
	    edge == RIGHT_EDGE || edge == FACES);
  CH_assert(guardcell == 0 || guardcell == 1);

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


box_info_t Chombo_Uniform_Grid::GetBoxInfo(const int blkID,
					   const int gridStruct)
  const
{
  box_info_t box_info;

  // Definition of FLASH UG
  CH_assert (blkID == 0);
  CH_assert (m_ugBoxLayout.numBoxes(procID()) == 1);
  CH_assert (gridStruct == CENTER || gridStruct == SCRATCH_CTR); //DEV TEMP!!!

  // Zero all data for the case where NDIM < MDIM
  for (int i=0; i<MDIM; ++i)
  {
    box_info.deltas[i] = 0.0;
    box_info.bsize[i] = 0.0;
    box_info.coord[i] = 0.0;
    box_info.lowBndBox[i] = 0.0;
    box_info.highBndBox[i] = 0.0;
    box_info.lowLimits[i] = 0;
    box_info.highLimits[i] = 0;
    box_info.guardcells[i] = 0;
    box_info.corner[i] = 0;
    box_info.stride[i] = 0;
    box_info.isNextToLowDomain[i] = FLASH_TRUE;
    box_info.isNextToHighDomain[i] = FLASH_TRUE;
  }
  box_info.lrefine = 0;

  DataIterator dit = m_ugBoxLayout.dataIterator();
  const Box& b = m_ugBoxLayout[dit()];
  const IntVect& small = b.smallEnd(); // e.g. (32,0)
  const IntVect& big = b.bigEnd(); // e.g. (63,31)
  const IntVect& size = b.size(); // e.g. (32,32)

  // I want to access the ProblemDomain index information.
  // This requires a box and its smallEnd and BigEnd methods.
  const Box& domainBox = m_computationalDomain.domainBox();
  const IntVect& small_domain = domainBox.smallEnd();
  const IntVect& big_domain = domainBox.bigEnd();

  for (int i=0; i<SpaceDim; ++i)
  {
    box_info.deltas[i] = m_deltas[i];
    box_info.bsize[i] = size[i] * m_deltas[i];
    box_info.lowBndBox[i] = small[i] * m_deltas[i];
    box_info.highBndBox[i] = (big[i] + 1) * m_deltas[i];
    box_info.coord[i] = box_info.lowBndBox[i] +
      ((box_info.highBndBox[i] - box_info.lowBndBox[i]) / 2.0);
    box_info.lowLimits[i] = small[i];
    box_info.highLimits[i] = big[i];
    box_info.guardcells[i] = m_ghostVect[i];
    box_info.corner[i] = small[i];
    box_info.stride[i] = 1;
    box_info.isNextToLowDomain[i] =
      (small[i] == small_domain[i] ? FLASH_TRUE : FLASH_FALSE);
    box_info.isNextToHighDomain[i] =
      (big[i] == big_domain[i] ? FLASH_TRUE : FLASH_FALSE);
  }

  return box_info;
}


void * Chombo_Uniform_Grid::GetDataPtr(const int blkID,
				       const int gridStruct)
{ 
  // Definition of FLASH UG
  CH_assert (blkID == 0);
  CH_assert (m_ugBoxLayout.numBoxes(procID()) == 1);

  LevelData<FArrayBox> aliasFab;
  switch (gridStruct) {
  case CENTER:
    CH_assert (m_ctrFab.isDefined());
    aliasLevelData(aliasFab, &m_ctrFab, m_ctrFab.interval());
    break;
  case SCRATCH_CTR:
    CH_assert (m_scratchCtrFab.isDefined());
    aliasLevelData(aliasFab, &m_scratchCtrFab, m_scratchCtrFab.interval());
    break;
  default: 
    MayDay::Error("[GetDataPtr]: Invalid mesh type");
  }

  DataIterator dit = aliasFab.dataIterator();
  FArrayBox& localBlkData = aliasFab[dit];
  return (void*) localBlkData.dataPtr(0);
}


// Fill ghost cells with valid data
void Chombo_Uniform_Grid::FillGuardCells()
{
  m_ctrFab.exchange();
}


std::vector<int> Chombo_Uniform_Grid::GetBlockIDs(const int blkType,
						  const int level)
  const
{
  // The indices of vBlk are used to access a cached box lookup table.
  // It gives FLASH fast access to Chombo boxes.  The indices give
  // access to all the boxes that match blkType (and level) criteria.
  const int gridStruct = CENTER;
  std::vector<int> vBlk;
  
  if (blkType == ALL_BLKS || blkType == LEAF)
  {
    vBlk.push_back(0);
  }
  else if (blkType == REFINEMENT)
  {
    CH_assert (level == 0);
    vBlk.push_back(0);
  }
  else
  {
    // Return all indices for boxes next to chosen boundaries
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

    box_info_t boxInfo = GetBoxInfo(0, gridStruct);

    // Test whether block i is touching dimension j boundary
    for (int j=s; j<=e; ++j)
    {
      if (boxInfo.isNextToLowDomain[j] == FLASH_TRUE ||
	  boxInfo.isNextToHighDomain[j] == FLASH_TRUE)
      {
	vBlk.push_back(0);
	break; // So we do not add the same block multiple times
      }
    }
  }
  return vBlk;
}

bool Chombo_Uniform_Grid::IsInitialRefinementDone()
  const
{
  return true;
}

void Chombo_Uniform_Grid::BuildInitialGrid()
{
}


void Chombo_Uniform_Grid::RefineInitialGrid()
{
}


void Chombo_Uniform_Grid::Regrid(const int baseLevel)
{
}


void Chombo_Uniform_Grid::AverageLevelData()
{
}

void Chombo_Uniform_Grid::ZeroFluxData()
{
}

void Chombo_Uniform_Grid::PutFluxData(const void *pFluxes,
				      const Real dt,
				      const int blkID,
				      const int axis)
{
}

void Chombo_Uniform_Grid::Reflux()
{
}

void Chombo_Uniform_Grid::PostTimeStep()
{
}

void Chombo_Uniform_Grid::FinalizeInitialGrid()
{
}


void Chombo_Uniform_Grid::PrintMeshNames()
  const 
{
  const int mesh_def[] =
    { CENTER, FACEX, FACEY, FACEZ, SCRATCH, SCRATCH_CTR,
      SCRATCH_FACEX, SCRATCH_FACEY, SCRATCH_FACEZ};
  const std::string mesh_str[] = {
    "ctr", "facex", "facey", "facez", "scratch", "scratchCtr",
    "scratchFacex", "scratchFacey", "scratchFacez"};
  const int num = sizeof(mesh_def) / sizeof(mesh_def[0]);

  // First print out the mesh strings for the keys we know.
  pout() << endl;
  for (int i=0; i<num; ++i)
  {
    const int def = mesh_def[i];
    const std::string str = mesh_str[i];
    mesh_info_t::const_iterator it = m_meshInfo.find(def);
    if (it != m_meshInfo.end())
    {
      pout() << "Mesh struct '" << str << "' corresponds to defined value " <<
	def << " and has " << it->second.size() << " elements " << std::endl;
      for (vector<string>::const_iterator itv = it->second.begin();
	   itv != it->second.end(); ++itv)
      {
	pout() << *itv << endl;
      }
    }
  }

  // Now print out all mesh strings in the container.
  pout() << endl;
  for (mesh_info_t::const_iterator it = m_meshInfo.begin(); 
       it != m_meshInfo.end(); ++it)
  {
    pout() << "Mesh defined value " << it->first << " has " << 
      it->second.size() << " elements " << endl;
    for (vector<string>::const_iterator itv = it->second.begin();
	 itv != it->second.end(); ++itv)
    {
      pout() << *itv << endl;
    }
  }
}
