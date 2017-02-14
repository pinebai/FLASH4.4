#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "mpi.h"
#include "flash_ptio.h"



/*I'll assume that MPI has already been setup.*/

/* readFlashParticles
*  inputs: filename, nprocs, mype
*  outputs: localParticles, localnp, tagIndex
*
*  Notes: make sure you deallocate the localParticles array!
*/         

int readFlashParticles(hid_t file_identifier, double* localParticles, 
                       int localnp, int offset, int npart_props){
    
    hid_t    dataspace, dataset, memspace;
    herr_t   status;

    int      rank, i;
    hsize_t  dimens_2d[2], start_2d[2], count_2d[2], stride_2d[2];
    
   
    /*read in chunk of particles*/

    dataset = H5Dopen(file_identifier, "tracer particles");
    if(dataset < 0){
	fprintf(stderr, "Error Opening Dataset!\n");
	return DATASET_OPEN_FAIL;
    }

    dataspace = H5Dget_space(dataset);
    if(dataspace < 0){
	fprintf(stderr, "Error Setting up Dataspace!");
	return DATASPACE_SELECT_FAIL;
    }


    start_2d[0] = (hsize_t) offset;
    start_2d[1] = (hsize_t) 0;

    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0] = localnp;
    count_2d[1] = npart_props; 

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, stride_2d,
				  count_2d, NULL);
    if (status < 0){
	fprintf(stderr,"Error Selecting Hyperslab!\n");
	return HYPERSLAB_SELECT_FAIL;
    }
    if(localnp > 0){

	rank = 2;
	dimens_2d[0] = localnp;
	dimens_2d[1] = npart_props;

	
       /* printf("MEMSPACE: %d %d %d \n", rank, localnp, npart_props);*/
        memspace = H5Screate_simple(rank, dimens_2d, NULL);
	if(memspace < 0){
	    fprintf(stderr,"Error selecting memspace!\n");
	    return MEMSPACE_SELECT_FAIL;
	}

	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
	                 H5P_DEFAULT, localParticles);
	if (status < 0){
	    printf("Error Reading In Particle Data!");
	    return DATA_READ_ERROR;
	}
	status = H5Sclose(memspace);
	if (status < 0){
	    printf("MEMSPACE\n");
	    exit(1);
	}
    }
   
    status = H5Sclose(dataspace);
   if (status < 0){
        printf("DATASPACE\n");
	exit(1);
    }
  
    status = H5Dclose(dataset);
    if (status < 0){
        printf("DATASET\n");
	exit(1);
    }
   
   return NORMAL_STATUS;
}


/*find local number of particles so you can pass in appropriately sized 
  localParticles array*/

int readFlashLocalNP(hid_t file_identifier, int *localnp, int *maxnp, 
		     int* start, int* end, int *tagIndex, int *npart_props,
                     int nprocs, int mype, MPI_Comm comm){

    hid_t dataspace, dataset, memspace, string_type;
    herr_t status;

    int rank, i;

    hsize_t start_1d, stride_1d, count_1d, dimens_1d;
    hsize_t start_2d[2], stride_2d[2], count_2d[2], dimens_2d[2];
    

    int ierr;
    int globalnp;
    int *np;

    hsize_t  maximum_dims[10];
    hsize_t  dataspace_dims[10];
    /*find the number of particle properties and which one is the tag*/
    char part_names[80][25];

    if(MASTERPE == mype){
	
	/*open the dataset*/
	/* read in whole thing, sum and broadcast to everyone */
	dataset = H5Dopen(file_identifier, "localnp");
	if(dataset < 0){
	    return DATASET_OPEN_FAIL;
	}
    
	dataspace = H5Dget_space(dataset);
	if(dataspace < 0){
	   return DATASPACE_SELECT_FAIL;
	}
	start_1d = (hsize_t)0;
	
	stride_1d = (hsize_t)1;

	count_1d = H5Sget_simple_extent_npoints(dataspace);

	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
				    &stride_1d, &count_1d, NULL);
       
	if(status < 0){
	  return HYPERSLAB_SELECT_FAIL;
	}

	rank = 1;
	dimens_1d = count_1d;

	memspace = H5Screate_simple(rank, &dimens_1d, NULL);
	if(memspace < 0){
	  return MEMSPACE_SELECT_FAIL;
	}

	np = (int *)malloc(sizeof(int)*count_1d);
	
	status  = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
			H5P_DEFAULT, np);
        
        globalnp = 0;
        for(i = 0; i < count_1d; ++i){
	   globalnp += np[i];
	}
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
    }

    MPI_Bcast(&count_1d, 1, MPI_INTEGER, MASTERPE, comm);
    MPI_Bcast(&globalnp, 1, MPI_INTEGER, MASTERPE, comm);

    /*everyone figure out how many particles we should have, MPE_Decomp1d 
      does exactly this*/

    MPE_Decomp1d(globalnp, nprocs, mype, start, end);
     /*printf("MyPE: %d: decomp1d args: %d %d %d %d %d\n", mype,
	    globalnp, nprocs, mype, *start, *end);   */
    *localnp = *end - *start + 1;
    *start -= 1; /*adjust for C's zero indexing */
    

    /*We also have to get particle metadata (number of particles and the 
      location of the tag index*/



    dataset = H5Dopen(file_identifier, "particle names");
    if(dataset < 0){
        return DATASET_OPEN_FAIL;
    }
    
    dataspace = H5Dget_space(dataset);
    if(dataspace < 0){
        return DATASPACE_SELECT_FAIL;
    }
    
    start_2d[0] = (hsize_t)0;
    start_2d[1] = (hsize_t)0;
    
    stride_2d[0] = (hsize_t)1;
    stride_2d[1] = (hsize_t)1;
    
    H5Sget_simple_extent_dims(dataspace,dataspace_dims, maximum_dims);
    *npart_props = count_2d[0] = dataspace_dims[0];
    count_2d[1] = 1;

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
                                 stride_2d, count_2d, NULL);
    
    if(status < 0){
        fprintf(stderr, "Hyperslab Selection Failed: Particle names!\n");
        return HYPERSLAB_SELECT_FAIL;
    }
    
    rank = 2;
    dimens_2d[0] = *npart_props;
    dimens_2d[1] = 1;
    
    /*part_names = (char**)malloc(sizeof(char*)*(*npart_props));
    for(i = 0; i < *npart_props; ++i)
        part_names[i] = (char *)malloc(sizeof(char)*25);*/

    /*must manually set the string size for now*/
    string_type = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(string_type, 25);

    memspace = H5Screate_simple(rank, dimens_2d, NULL);
    if(memspace < 0){
	fprintf(stderr, "Memspace Creation for particle names failed!\n");
	return MEMSPACE_SELECT_FAIL;
    }
    
    status = H5Dread(dataset, string_type, memspace, dataspace, H5P_DEFAULT,
                     part_names);
    if(status < 0){
      fprintf(stderr, "Error Reading particle name data!");
      return DATA_READ_ERROR;
      }

    for(i = 0; i < *npart_props; ++i){
      if(!(strncmp(part_names[i], "tag", 3))){
        *tagIndex = i;
        /*break;*/ 
      }
    }
 
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    MPI_Barrier(MPI_COMM_WORLD);
    return NORMAL_STATUS;
}

int flashFileOpen(char* filename, hid_t* file_identifier){

    int ierr;
    MPI_Info fileInfoTemplate;
    hid_t acc_template;

    /*set up file access*/
    ierr = MPI_Info_create(&fileInfoTemplate);

    acc_template = H5Pcreate(H5P_FILE_ACCESS);
    /*open for a parallel read*/

    ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);

    *file_identifier = H5Fopen(filename, H5F_ACC_RDWR, acc_template);

    ierr = H5Pclose(acc_template);
    
    if(*file_identifier < 0){
        fprintf(stderr, "Error opening file, %s!\nError returned is %d!",
	         filename, *file_identifier);
	return FILE_OPEN_FAIL;
    }
    
    return NORMAL_STATUS;
}

void flashFileClose(hid_t* file_identifier){

    MPI_Barrier(MPI_COMM_WORLD);
    H5Fclose(*file_identifier);
    
    return;
}

		     
