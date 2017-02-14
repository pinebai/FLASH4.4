#include "DataBufferComparer.hpp"

DataBufferComparer::DataBufferComparer(const MPI_Comm &comm)
{
  MPI_Comm_rank(comm, &m_myPE);
  MPI_Comm_size(comm, &m_numProcs);
  m_comm = comm; 
  m_comparisonFailed = false;
}


DataBufferComparer::~DataBufferComparer()
{
  return;
}


double DataBufferComparer::max(const double &a, const double &b) const {
  return a > b ? a : b;
}


double DataBufferComparer::norm(const double &a, const double &b) const {
  return fabs(2*(a-b))/(fabs(a+b) > 1e-99 ? fabs(a+b) : 1e-99); 
}


void DataBufferComparer::DoComparison(const std::string &varName, 
				      const double * const buf1, 
				      const double * const buf2, 
				      const size_t numDataElements) const
{
  const size_t TOT_IDX = 3;
  const size_t MIN_IDX = 0, MAX_IDX = 1, ABS_IDX = 2;

  struct {
    double error;
    int pe;
  } localStat[TOT_IDX], globalStat[TOT_IDX];
  double localCellTuple[TOT_IDX*2], globalCellTuple[TOT_IDX*2];

  double error = 0.0, globalError;
  long long globalNumMismatch = 0, globalNumDataElements = 0, 
    localNumMismatch = 0, localNumDataElements = 0;
  size_t numMismatch = 0;
  int flagLocalFailure = 0, flagGlobalFailure = 0;
  bool hasInitialValue[TOT_IDX];
  double min1, min2, max1, max2;

  /* DBL_MAX and DBL_MIN are defined in cfloat */
  localStat[MAX_IDX].error = 0.0;     localStat[MAX_IDX].pe = m_myPE;
  localStat[MIN_IDX].error = DBL_MAX; localStat[MIN_IDX].pe = m_myPE;
  localStat[ABS_IDX].error = 0.0;     localStat[ABS_IDX].pe = m_myPE;
  min1 = DBL_MAX;
  min2 = DBL_MAX;
  max1 = DBL_MIN;
  max2 = DBL_MIN;
  error = DBL_MIN;

  /* Compare the files' data */
  if ((buf1 != NULL) && (buf2 != NULL)) {
    for(size_t i=0; i<numDataElements; ++i) {
      if (buf1[i] != buf2[i]) {
	error = norm(buf1[i], buf2[i]);

	if (error > localStat[MAX_IDX].error) {
	  localStat[MAX_IDX].error = error;
	  localCellTuple[MAX_IDX] = buf1[i];
	  localCellTuple[MAX_IDX+1] = buf2[i];
	}

	if (error < localStat[MIN_IDX].error) {
	  localStat[MIN_IDX].error = error;
	  localCellTuple[MIN_IDX] = buf1[i];
	  localCellTuple[MIN_IDX+1] = buf2[i];
	}

	error = fabs(buf1[i] - buf2[i]);
	if (error > localStat[ABS_IDX].error){
	  localStat[ABS_IDX].error = error;
	  localCellTuple[ABS_IDX] = buf1[i];
	  localCellTuple[ABS_IDX+1] = buf2[i];
	}

	numMismatch = numMismatch + 1;
	flagLocalFailure = 1;
      }
      max1 = buf1[i] > max1 ? buf1[i] : max1;
      max2 = buf2[i] > max2 ? buf2[i] : max2;
      min1 = buf1[i] < min1 ? buf1[i] : min1;
      min2 = buf2[i] < min2 ? buf2[i] : min2;
    }

    error = fabs(max1);
    error = fabs(min1) > error ? fabs(min1) : error;
    error = fabs(max2) > error ? fabs(max2) : error;
    error = fabs(min2) > error ? fabs(min2) : error;
    error = error > 1.e-99 ? error : 1.e-99;
  }


#ifdef DEBUG
  std::fstream File;
  std::stringstream ss;
  ss << std::setw(4) << std::setfill('0') << m_myPE;
  std::string localLogFile = "debug." + ss.str();

  File.open(localLogFile.c_str(), std::ios::out | std::ios::app);
  File << varName << " - has " << numMismatch << " cell mismatches from " <<
    numDataElements << " cells" << std::endl;
  File.close();
#endif


  MPI_Allreduce(&flagLocalFailure, &flagGlobalFailure, 1, MPI_INT, MPI_LOR, m_comm);

  if (flagGlobalFailure) {

    m_comparisonFailed = true;

    /* Make sure all integers that are used in the ALL_REDUCE operations are of 
       type long long.  I had problems when I mixed integer types on BG/P.  Also
       performing an ALL_REDUCE on type unsigned long long failed on BG/P, but
       this is probably an MPI bug.*/
    localNumMismatch = (long long) numMismatch;
    localNumDataElements = (long long) numDataElements;

    MPI_Allreduce(&localNumMismatch, &globalNumMismatch, 1, 
		  MPI_LONG_LONG, MPI_SUM, m_comm);
    MPI_Allreduce(&localNumDataElements, &globalNumDataElements, 1, 
		  MPI_LONG_LONG, MPI_SUM, m_comm);

    MPI_Allreduce(&localStat[MAX_IDX], &globalStat[MAX_IDX], 1,
		  MPI_DOUBLE_INT, MPI_MAXLOC, m_comm);
    MPI_Allreduce(&localStat[MIN_IDX], &globalStat[MIN_IDX], 1,
		  MPI_DOUBLE_INT, MPI_MINLOC, m_comm);
    MPI_Allreduce(&localStat[ABS_IDX], &globalStat[ABS_IDX], 1,
		  MPI_DOUBLE_INT, MPI_MAXLOC, m_comm);

    /* The processor which has the most significant error now broadcasts the
       corresponding cell tuple */
    for (size_t i=0; i<TOT_IDX; ++i) {
      size_t idx = i*2;
      if (globalStat[i].pe == localStat[i].pe) {
	globalCellTuple[idx] = localCellTuple[idx];
	globalCellTuple[idx+1] = localCellTuple[idx+1];
      }
      MPI_Bcast(&globalCellTuple[idx], 2, MPI_DOUBLE, globalStat[i].pe, m_comm);
    }

    MPI_Allreduce(&error, &globalError, 1, MPI_DOUBLE, MPI_MAX, m_comm);

    if (m_myPE == 0) {
      std::cout << std::setprecision(4) << "!!!!!! FILE MISMATCH (" << varName << 
	") !!!!!!\nNumber of cell mismatches " << globalNumMismatch << 
	" out of " << globalNumDataElements << 
	"\nMax error " << globalStat[MAX_IDX].error << 
	" (file1 " << globalCellTuple[MAX_IDX] <<
	", file2 " << globalCellTuple[MAX_IDX+1] << ")" <<
	"\nMin error " << globalStat[MIN_IDX].error << 
	" (file1 " << globalCellTuple[MIN_IDX] <<
	", file2 " << globalCellTuple[MIN_IDX+1] << ")" <<
	"\nAbs error " << globalStat[ABS_IDX].error << 
	" (file1 " << globalCellTuple[ABS_IDX] <<
	", file2 " << globalCellTuple[ABS_IDX+1] << ")" <<
	"\nMag error " << globalStat[ABS_IDX].error / globalError <<
	std::endl;
    }
  } else {
    if (m_myPE == 0) {
      std::cout << "All elements match for " << varName << std::endl;
    }
  }
}


bool DataBufferComparer::didComparisonFail(){
  return m_comparisonFailed;
}
