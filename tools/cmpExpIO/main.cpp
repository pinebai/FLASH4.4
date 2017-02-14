#include "mpi.h"
#include <algorithm>
#include "FlashFileFactory.hpp"
#include "DataBufferComparer.hpp"


int main(int argc, char *argv[])
{
  std::vector<std::string> fileStrings, variableNames, v1, v2;
  FlashFile *file1 = NULL, *file2 = NULL;
  const double *buf1 = NULL, *buf2 = NULL;
  size_t n1, n2;
  int myPE, numProcs, i;
  std::string infoStr;
  bool useCollectiveIO = true;

  // The option handling is ugly.  I would use getopt, but I've
  // already seen sfocu crash on BG/Q in getopt code.
  if (argc == 3 || argc == 4) {
    for (i=1; i<argc; ++i) {
      if (std::string(argv[i]) == "-i") {
	if (!useCollectiveIO) {
	  std::cerr << "Multiple -i options passed\n";
	  exit(EXIT_FAILURE);
	}
	useCollectiveIO = false;
      } else {
	fileStrings.push_back(std::string(argv[i]));
      }
    }
  } else {
    std::cerr << std::endl << "Usage: " << argv[0] << 
      " checkpoint_file1 checkpoint_file2 [-i]" << std::endl;
    exit(EXIT_FAILURE);
  }


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (myPE == 0) {
    std::string xferStr;
    if (useCollectiveIO) {
      xferStr = "collective parallel I/O";
    } else {
      xferStr = "independent parallel I/O";
    }
    std::cout << "Comparison tool running on " << numProcs << 
      " processors - using " << xferStr << "\n" <<
      argv[0] << " " << fileStrings[0] << " " << fileStrings[1] << std::endl;
  }
  FlashFileFactory flashFileFactory(myPE, numProcs, useCollectiveIO);


  /* Create collective transfer property list for parallel dataset 
     write and initialise HDF5 derived data types. */
  file1 = flashFileFactory.GetFlashFileInstance(fileStrings[0]);
  file2 = flashFileFactory.GetFlashFileInstance(fileStrings[1]);
  if ( (NULL == file1) || (NULL == file2) ) {
    std::cerr << "NULL factory file object returned.\n";
    exit(EXIT_FAILURE);
  }


  /* Check that the number of grid data elements are the same. */
  n1 = file1->GetNumberDataElements();
  n2 = file2->GetNumberDataElements();
  if (n1 != n2) {
    std::cerr <<  "File1:" << n1 << " File2:" << n2
	      << "\n!!!!!! FILE MISMATCH (number of grid points) !!!!!!\n";
    exit(EXIT_FAILURE);
  }


  /* Constuct a vector of unique variable names. */
  v1 = file1->GetAllVariableNames();
  v2 = file2->GetAllVariableNames();
  
  variableNames.reserve(v1.size() + v2.size());
  variableNames.insert(variableNames.end(), v1.begin(), v1.end());
  variableNames.insert(variableNames.end(), v2.begin(), v2.end());

  std::sort(variableNames.begin(), variableNames.end());
  std::vector<std::string>::const_iterator endLocation = 
    std::unique(variableNames.begin(), variableNames.end());



  /* Loop over the unique string names only.  Extract data from files 
     for each variable and then compare. */
  DataBufferComparer dataBufferComparer(MPI_COMM_WORLD);

  for(std::vector<std::string>::const_iterator it = variableNames.begin(); 
      it != endLocation; ++it) {

    if (myPE == 0) {
      std::cout << "\nTesting " << *it << std::endl;
    }
    buf1 = file1->GetVariableFromFile(*it);
    buf2 = file2->GetVariableFromFile(*it);

    /* Note: We pass "n1" data elements since we know the number of elements 
       are the same in both files */
    dataBufferComparer.DoComparison(*it, buf1, buf2, n1);
  }

  if (myPE == 0) {
    if(dataBufferComparer.didComparisonFail())
      std::cout << "comparison of benchmark files yielded: FAILURE" << std::endl;
    else
      std::cout << "comparison of benchmark files yielded: SUCCESS" << std::endl;
  }

  delete file1;
  delete file2;

  MPI_Finalize();
  return 0;
}
