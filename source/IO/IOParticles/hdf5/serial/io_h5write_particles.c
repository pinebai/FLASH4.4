#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "constants.h"


int Driver_abortFlashC(char* message);


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5write_particles)(int* myPE,
                        hid_t* file_identifier,
                        int* totalparticles,
                        int* localnp,
                        int* npart_props,
                        int* particle_offset,
                        double* particles,
                        char particle_labels[][OUTPUT_PROP_LENGTH + 1])      /* particle labels */
{
  
  hid_t    dataspace, dataset, memspace, dataset_plist;
  herr_t   status, ierr;
  
  int      rank, i, j, string_size;
  hsize_t  dimens_2d[2], maxdims_1d;
  
  hsize_t start_2d[2];
  hsize_t  stride_2d[2], count_2d[2];
  
  /*the particle data type*/
  hid_t    particle_tid, string_type;
  
  char     labelstring[OUTPUT_PROP_LENGTH+1];


    
  if(*totalparticles > 0) {
    
    /* first write out the particle attribute names */
    
    rank = 2;
    dimens_2d[0] = (hsize_t) (*npart_props);
    dimens_2d[1] = 1;
    
    /*manually set the string size */
    string_size = OUTPUT_PROP_LENGTH;
    
    /*setup the datatype for this string length */
    string_type = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(string_type, string_size);
    if(status < 0){
      Driver_abortFlashC("Error: string\n");
    }
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    if(dataspace < 0){
      Driver_abortFlashC("Error: H5Screate_simple io_h5write_particles\n");
    }
    
    
    
    /* only need 1 proc to write particle labels one time */
    /* this if statement handles both the parallel and serial case */
    if((*particle_offset == 0) && (*localnp !=0)) {
      dataset   = H5Dcreate(*file_identifier, "particle names", 
			    string_type, dataspace, H5P_DEFAULT);
      if(dataset < 0){
	Driver_abortFlashC("Error: H5Dcreate particle names, io_h5write_particles\n");
      }     
      
      
      status    = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, 
			   H5P_DEFAULT, particle_labels);
      if(status < 0){
	Driver_abortFlashC("Error io_h5write_particle: H5Dwrite particle labels\n");
      }
      H5Dclose(dataset);  
    }
   
    H5Tclose(string_type);
    H5Sclose(dataspace);
    
  }
  

    
    
  /* then write out the particle data */  
  
  /* create the particle data space */
  rank       = 2;
  dimens_2d[0] = (hsize_t) (*totalparticles);
  dimens_2d[1] = (hsize_t) (*npart_props);
  
  
  
  dataspace  = H5Screate_simple(rank, dimens_2d, NULL);
  if(dataspace < 0) {
    Driver_abortFlashC("Error io_h5write_particles: H5Screate_simple tracer particles\n");
  }
  

  if(*localnp != 0){
  
    /* if the dataset has already been created then just open it */
    if ((*particle_offset != 0)) {
      dataset = H5Dopen(*file_identifier, "tracer particles"); 
      if(dataset < 0) {
	Driver_abortFlashC("Error: H5Dopen io_h5write_particles\n");
      }
      
    }else {
      /* create the dataset */ 
      dataset = H5Dcreate(*file_identifier, "tracer particles", H5T_NATIVE_DOUBLE,
			  dataspace, H5P_DEFAULT);
      if(dataset < 0) {
	Driver_abortFlashC("Error io_h5write_particles: H5Dcreate tracer particles\n");
      }
      
    }
  
    
    start_2d[0] = (hsize_t) (*particle_offset);
    start_2d[1] = 0;
    
    stride_2d[0] = (hsize_t)1;
    stride_2d[1] = (hsize_t)1;
    
    count_2d[0] = (hsize_t) (*localnp);
    count_2d[1] = (hsize_t) (*npart_props);
    
    
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
    stride_2d, count_2d, NULL);
    
    if(ierr < 0){
      printf("%s\n", "Error: unable to select hyperslab for particles dataspace");
      Driver_abortFlashC("Error: unable to select hyperslab for particles dataspace");
    }
    
    
    /* create the memory space */
    rank       = 2;
    dimens_2d[0] = (hsize_t) (*localnp);
    dimens_2d[1] = (hsize_t) (*npart_props);
    
    memspace   = H5Screate_simple(rank, dimens_2d, NULL);
    
    start_2d[0] = (hsize_t) 0;
    start_2d[1] = 0;
    
    stride_2d[0] = (hsize_t)1;
    stride_2d[1] = (hsize_t)1;
    
    count_2d[0] = (hsize_t) (*localnp);
    count_2d[1] = (hsize_t) (*npart_props);
    
    
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_2d,
				 stride_2d, count_2d, NULL); 
    
    if(status < 0) {
      printf("%s\n","Unable to write tracer particles");
      Driver_abortFlashC("Unable to write tracer particles");
    }
    
    
      
    /* write data to the dataset */ 
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
		      particles);
      
    
    if(status < 0) {
      printf("%s\n","Unable to write tracer particles");
      Driver_abortFlashC("Unable to write tracer particles");
    }
    /* release resources */
    status = H5Sclose(memspace);
    status = H5Dclose(dataset); 

  }
  
  
  /* release resources */
  status = H5Sclose(dataspace);

}
    


  








