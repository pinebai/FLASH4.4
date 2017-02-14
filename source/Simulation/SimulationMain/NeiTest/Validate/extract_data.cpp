#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "preflash_datafile.hpp"
#include "preflash_dataset.hpp"
#include "preflash_siminfo.hpp"
#include "preflash_utils_text.hpp"
#include "preflash_utils_file.hpp"
#include "preflash_utils_commandline.hpp"

using std::cout; using std::cerr; using std::endl;
using std::ofstream; using std::scientific;
using std::vector; using std::string; using std::list;

const std::string CreatePlottableLabel(const std::string &);

/* The purpose of this program is to extract data for every 
single variable in a given HDF5 input file.  The HDF5 file is assumed 
to have come from a 1D simulation.  We will produce a separate 
ASCII text file for each variable.  The output files will then 
be post-processed and graphed in a python script. */

int main(int argc, char * argv[])
{
  vector<string> argument_list;
  list<PreFlash::Utils::Option> option_list;

  PreFlash::Utils::parse_command_line(argc, argv, argument_list, option_list);

  const uint num_args = argument_list.size();
  const uint required_args = 1;

  if (num_args != required_args)
  {
    const string full_exec = argv[0];
    string exec_path, exec_name;

    PreFlash::Utils::split_name_path_file(full_exec, exec_path, exec_name);
    cerr << endl << "Usage: " << exec_name << " dataset_name" << endl;
    return EXIT_FAILURE;
  }
  

  /* We are not interested in any special options, we just want a single 
     argument containing the HDF5 file name */
  const string datafile_name = argument_list[0];
  cout << "We are working with the HDF5 file: " << datafile_name << endl;


  /* Read the HDF5 file */
  const uint buffer_size = 0;   /* Use default buffer size */
  const bool skip_particles = true;
  const PreFlash::File::DataFile dfile(datafile_name, buffer_size, skip_particles);

  /* NOTE: This is returning the value in an enum:
     FLASH2: 0, FLASH3: 1, FLASH3.1: 2 */
  const int flashVersion = dfile.get_version();
  cout << "Flash version is: " << flashVersion
       << " where F2=0, F3=1, F3.1=2." << endl;

  /* Print out various simulation information */
  const PreFlash::File::SimInfo & siminfo = dfile.get_sim_info();
  const uint dims = siminfo.get_dims();
  cout << "Dimensions: " << dims << endl;
  if (dims != 1) 
  {
    cerr << "This is a limited program only designed to work for 1D datasets" << endl;
    return EXIT_FAILURE;
  }


  const vector<double> & volume_min_bounds = siminfo.get_volume_minbounds();
  const vector<double> & volume_max_bounds = siminfo.get_volume_maxbounds();
  cout << "Volume min bounds: " << volume_min_bounds << endl;
  cout << "Volume max bounds: " << volume_max_bounds << endl;

  const double sim_time = siminfo.get_sim_time();
  cout << "Simulation time: " << sim_time << endl;

  const uint num_blocks = siminfo.get_num_blocks();
  cout << "Number of blocks: " << num_blocks << endl;
  cout << "Leaf blocks: " << siminfo.get_num_leaf_blocks() << endl;
  cout << "Base (top level) blocks: " << siminfo.get_base_block_dims() << endl;
  cout << "Number of cells in a block: " << siminfo.get_num_block_cells() << endl;


  /* Iterate over each leaf block */
  if (num_blocks > 0)
  {
    uint min_leaf_refine_level = 999;
    uint max_leaf_refine_level = 0;

    for (uint index = 0; index < num_blocks; index++)
    {
      if (siminfo.is_leaf(index))
      {
	const uint refine_level = siminfo.get_refine_level(index);

	if (refine_level < min_leaf_refine_level)
	  min_leaf_refine_level = refine_level;
	else if (refine_level > max_leaf_refine_level)
	  max_leaf_refine_level = refine_level;
      }
    }
    
    cout << endl << "Leaf block refinement level varies from " << min_leaf_refine_level << 
      " to " << max_leaf_refine_level << endl;
  }



  /* Iterate over each dataset, and write the cell values
     to unique files */
  const vector<string> & dataset_names = siminfo.get_mesh_variable_names();
  cout << endl << "Mesh variable names: " << dataset_names  << endl;


  const uint num_datasets = dataset_names.size();
  for (uint dataset_index = 0; dataset_index < num_datasets;
       dataset_index++)
  {

    const string & dset_name = dataset_names[dataset_index];
    cout << "Extracting data for variable: " << dset_name << endl;


    /* Create a datafile name for this unk variable.  Note that 
       we modify names so that they are the same as FLASH2 */
    std::string outFile;
    if (flashVersion > 0) {
      const std::string F2_dest_name = CreatePlottableLabel(dset_name);
      outFile = F2_dest_name + ".dat";
      cout << "FLASH2 compatibility... label this variable's data file: " << outFile << endl;
    } else {
      outFile = dset_name + ".dat";
    }
    ofstream outputStream(outFile.c_str());


    string padded_dset_name;   
    PreFlash::File::pad_variable_name(dset_name, padded_dset_name);
    const PreFlash::File::Dataset & dset = dfile.get_dataset(padded_dset_name);
    PreFlash::Block::BlockData<double> stored_data;


    /* Loop over each block for a given dataset */
    for (uint block_index = 0; block_index < num_blocks; block_index++)
    {
      if (siminfo.is_leaf(block_index))
      {
	dset.get_block_data(block_index, stored_data);

	const PreFlash::Block::BlockInfo & block_info 
	  = siminfo.get_block_info(block_index);

	const uint num_cells = block_info.get_num_cells();
	vector<double> cell_min_coords(dims);
	vector<double> cell_max_coords(dims);
	vector<double> cell_given_position(dims);


	/* Loop over each cell for a given block */	
	for (uint cell_index = 0; cell_index < num_cells; cell_index++)
	{
	  const double cell_value = stored_data[cell_index];
	  block_info.get_cell_bounds(cell_index, cell_min_coords, 
				     cell_max_coords);
	  
	  for (uint idim = 0; idim < dims; ++idim) 
	    cell_given_position[idim] = 
	      (cell_min_coords[idim] + cell_max_coords[idim]) / 2.0;
	  
	  outputStream << scientific << cell_given_position << " " << cell_value << endl;
	}
      }
    }
    outputStream.close();
    cout << "Finished writing data for variable: " << dset_name << endl << endl;
  }
  
  return EXIT_SUCCESS;
}


//Create labels like we use in FLASH2.
const std::string CreatePlottableLabel(const std::string & dset_name)
{

  std::string numericString="", characterString="";
  int numericValue;


  for (std::string::const_iterator i = dset_name.begin(); 
       i != dset_name.end(); ++i) 
  {
    if (isdigit(*i)) {
      numericString += *i;
    } else {
      characterString += *i;
    }
  }


  if (numericString.empty()) {
    numericString = "01";
  } else {

    //Convert numericString to an integer.  Add one. 
    std::istringstream ssInput(numericString);
    ssInput >> numericValue;
    numericValue += 1;
    
    //Copy new value back to original string.
    std::ostringstream ssOutput;
    ssOutput.width(2);
    ssOutput.fill('0');
    ssOutput << numericValue;
    numericString = ssOutput.str();    
  }

  if (characterString.size() == 1) {
    characterString += " ";
  }

  const std::string newString = characterString + numericString;
  return newString;
}
