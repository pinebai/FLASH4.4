#include <string.h>
#include <stdlib.h>
#include "hdf5.h"
#include "flash_reader.h"
#include "options.h"
#include "namecmp.h"

int FR_get_data_names_HDF5(hid_t filename, char names[FR_MAXVARS][FR_VAR_STRING_SIZE+1], int *var_type, int *num_vars, int *ncells_vec);

/* if mem_type is 0, use the one returned by H5Dget_type */
int FR_slurp_HDF5(hid_t handle, hid_t mem_type, char *name, int dim, void *target){
  hid_t dataset, dataset_type, dataspace, memspace;
  herr_t status;
  int auto_type = 0;
 
  /*for special handling*/
  hsize_t dimens_3d[3],maxdimens_3d[3], mem_dimens_3d[3];
  hsize_t count_3d[3], start_3d[3], stride_3d[3];
  
  hsize_t dimens_2d[2], maxdimens_2d[2], mem_dimens_2d[2];
  hsize_t count_2d[2], start_2d[2], stride_2d[2];
  int rank;
  


  if (!mem_type)
    auto_type = 1;

  dataset = H5Dopen(handle, name);
  if (dataset < 0) {
    printf("Error reading \"%s\": dataset not found\n", name);
    return 1; 
  }

  /* Not sure how safe it is to get mem_type from the dataset, so
   I'm trying this for a while.
   Eventually eliminate the mem_type argument.
   */
  if (auto_type) { /* get type from dataset */
    mem_type = H5Dget_type(dataset);
    dataset_type = mem_type;
  } else { /*just check type */
    dataset_type = H5Dget_type(dataset);
    if (H5Tget_class(mem_type) != H5Tget_class(dataset_type))
      printf("\"%s\" has an unexpected data type\n", name);
  }
  
  /*This has to change for bounding box*/
  if(strncmp("bounding box", name, 12) == 0){
    /*make sure memory space is as though this is 3d*/
    dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, dimens_3d, maxdimens_3d);
    
    
    /*we have to make sure we are going to correctly map to memory*/
    
    start_3d[0] = 0;
    start_3d[1] = 0;
    start_3d[2] = 0;
    
    stride_3d[0] = 1;
    stride_3d[1] = 1;
    stride_3d[2] = 1;
      
    count_3d[0] = dimens_3d[0];
    count_3d[1] = (hsize_t)dim;
    count_3d[2] = dimens_3d[2];
    
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
                        start_3d, stride_3d, count_3d, NULL);
    
    rank = 3;
    mem_dimens_3d[0] = dimens_3d[0];
    mem_dimens_3d[1] = (hsize_t)dim;
    mem_dimens_3d[2] = dimens_3d[2];
    
    memspace = H5Screate_simple(rank, mem_dimens_3d, NULL);
    
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                        start_3d, stride_3d, count_3d, NULL);
    
    status = H5Dread(dataset, mem_type, memspace, dataspace, H5P_DEFAULT, target);
    if (status < 0) {
      printf("Error reading \"%s\"\n", name);
      return 1; 
    }  
    H5Sclose(memspace);
    H5Sclose(dataspace);
  
  }/*end if bounding box*/
  
  else if ((strncmp(name, "coordinates", 11) == 0) || (strncmp("block size", name, 10) == 0)){
    dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
    
    
    /*we have to make sure we are going to correctly map to memory*/
    
    start_2d[0] = 0;
    start_2d[1] = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0] = dimens_2d[0];
    count_2d[1] = (hsize_t)dim;
    
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
                        start_2d, stride_2d, count_2d, NULL);
    
    rank = 2;
    mem_dimens_2d[0] = dimens_2d[0];
    mem_dimens_2d[1] = (hsize_t)dim;
        
    memspace = H5Screate_simple(rank, mem_dimens_2d, NULL);
    
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                        start_2d, stride_2d, count_2d, NULL);
    
    status = H5Dread(dataset, mem_type, memspace, dataspace, H5P_DEFAULT, target);
    if (status < 0) {
      printf("Error reading \"%s\"\n", name);
      return 1; 
    }  
    H5Sclose(memspace);
    H5Sclose(dataspace);
    
  } /*end coordinates handling*/
  
  
  /*normal slurp*/
  else{
    status = H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, target);
    if (status < 0) {
      printf("Error reading \"%s\"\n", name);
      return 1; 
    }
  }
  H5Tclose(dataset_type);
  H5Dclose(dataset);
  return 0; 
}

FR_ParticlesAllProps *FR_GetParticlesAllProps_HDF5(FR_File *file, int startingAt, 
						   int atMost, int *numberGotten) {
  FR_ParticlesAllProps *particles;
  int i;
  herr_t   status;
  hsize_t  dimens_1d, maxdims_1d;
  hssize_t start_2d[2];
  hsize_t  stride_2d[2], count_2d[2];
  hsize_t dims_out[2];
  int rank;
  hid_t dataset, dataspace, memspace;
  int to_read;
  hid_t particle_tid;
  
  
  if (file->totalparticles > 0) {
    if (atMost < (file->totalparticles - startingAt)) {
      to_read = atMost;
    } else {
      to_read = file->totalparticles - startingAt;
    }

    if (to_read > 0) {
      particles = (FR_ParticlesAllProps *) malloc(sizeof(FR_ParticlesAllProps));

      particles->numIntProps = file->numIntPartProps;
      particles->numRealProps = file->numRealPartProps;

      particles->realProps = (double *) malloc(sizeof(double) * to_read * particles->numRealProps);
      particles->intProps = (int *) malloc(sizeof(int) * to_read * particles->numIntProps);

      particles->realPropsNames = (char **) malloc(sizeof(char *) * particles->numRealProps);
      particles->intPropsNames  = (char **) malloc(sizeof(char *) * particles->numIntProps);

      
      for (i = 0; i < particles->numRealProps; i++) {
	particles->realPropsNames[i] = (char *) malloc(sizeof(char) * (FR_PART_PROP_STRING_SIZE+1));
	strcpy(particles->realPropsNames[i], file->realPartPropNames[i]);
      }
      

      for (i = 0; i < particles->numIntProps; i++) {
	particles->intPropsNames[i] = (char *) malloc(sizeof(char) * (FR_PART_PROP_STRING_SIZE+1));
	strcpy(particles->intPropsNames[i], file->intPartPropNames[i]);
      }

      dataset = H5Dopen(*((hid_t *)file->handle), "tracer particles");

      dataspace = H5Dget_space(dataset);

      rank       = 1;
      start_2d[0]   = startingAt;
      start_2d[1]   = 0;
      count_2d[0]   = (hsize_t) (to_read);
      count_2d[1]   = (hsize_t) (file->numRealPartProps);
      maxdims_1d = (hsize_t) (file->totalparticles * file->numRealPartProps);
      dimens_1d  = (hsize_t) (to_read*file->numRealPartProps);

      status     = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
                                       NULL, count_2d, NULL);
      
      memspace   = H5Screate_simple(rank, &dimens_1d, &maxdims_1d);

      /* read data from the dataset */
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
		       particles->realProps);

      if (status < 0) {
	printf("Error reading particles from file '%s'\n", file->filename);
      }
      else {
	(*numberGotten) = to_read;
      }

      H5Dclose(dataset);
      H5Sclose(dataspace);
      H5Sclose(memspace);
    
      *numberGotten = to_read;
      particles->numParticles = to_read;
    } 
  } else {
    *numberGotten = 0;
    particles = NULL;
  }
  return particles;
}

/*
 * GetNumParticles: look at the particles entry and get the size
 *                of the dataset based on the dimensions of the array
 *
 * input:  file handle          (integer)
 *
 * output: number of particles  (integer)
 *
 * return value: -1 indicates an error getting the dimension
 *
 */
void FR_GetNumParticles_HDF5(hid_t *file_identifier,
		      hid_t *numParticles)
{

  hid_t dataset, dataspace;
  hsize_t maximum_dims[10];
  hsize_t dataspace_dims[10];
  herr_t (*old_func)(void*);
  void *old_client_data;
  hid_t partTypeId;

  int numElements;

  /* temporarily turn off error handling, to probe for the dataset's existence */
  H5Eget_auto(&old_func, &old_client_data);
  H5Eset_auto(NULL, NULL);

  /* get the dataset ID for the tracer partlcies record, and the dataspace in
     this record */
  dataset = H5Dopen(*file_identifier, "tracer particles");

  /* restore the error handling */
  H5Eset_auto(old_func, old_client_data);

  if (dataset >= 0) {
    dataspace = H5Dget_space(dataset);

    /* get the dimensions of the dataspace */    
    H5Sget_simple_extent_dims(dataspace, dataspace_dims, maximum_dims);

    *numParticles = dataspace_dims[0];

    H5Dclose(dataset);
    H5Sclose(dataspace);

  } else {
    *numParticles = 0;
  }
}

/* from io/hdf5v1.4_serial_single_compound/hdf5_flash.h FIXME */ 
/* commented out in favor of int_scalars_t
typedef struct sim_params_t {
  int total_blocks;
  int nsteps;
  int nxb;
  int nyb;
  int nzb;
  double time; 
  double timestep;
  double redshift;
} sim_params_t;
*/
/*
typedef struct int_scalars_t {
  int nstep;
  int nbegin;
  int nxb;
  int nyb;
  int nzb;
  int globalnumblocks;
} int_scalars_t;
*/

typedef struct int_list_t {
  /* DEV - hard code for now the length of 'name' - NTT */
  char name[MAX_STRING_LENGTH];
  int value;
} int_list_t;

/* The strings in the 'name' field of the 'integer scalars' dataset
   have a constant length of 80 characters where the meaningful part
   of the word is at the beginning and the rest is padded out with
   spaces. Here we only compare the meaningful part of the word with
   the name we are looking for, ignoring the whitespace. */
int specialcmp(const char *a, const char *b) {
  for ( ; (*a) && (*a != ' '); a++, b++) {
    if ((! *b) || (*a != *b)) {
      return 0;
    }
  }
  /* if b has not terminated, it's not a match */
  if (*b) return 0;
  else return 1;
}


FR_File *FR_open_HDF5_Pmesh(hid_t handle,FR_File *out, options_t opts);
extern FR_File *FR_open_HDF5_Chombo(hid_t handle, FR_File* out, options_t opts);


FR_File *FR_open_HDF5(char *filename, options_t opts){
	FR_File *out;

	hid_t handle;
	hid_t gid;
	int status;

	/*Turn off automatic error printing*/
	status=H5Eset_auto(NULL, NULL);

  	/* Open file */
  	out = (FR_File *) malloc(sizeof(FR_File));
 	 out->handle = malloc(sizeof(hid_t));
 	 strcpy(out->filename, filename);

  	handle = H5Fopen(out->filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  	*((hid_t *)out->handle) = handle;



  	if (handle< 0)
    		return NULL;

  	/*Determine if this is a Chombo or Paramesh file */ 
	if( (gid= H5Gopen(handle, "Chombo_global")) < 0){
  		out->format = FR_HDF5_PMESH;
		out=FR_open_HDF5_Pmesh(handle, out, opts);
	}else{
		out->format = FR_HDF5_CHOMBO;
		out=FR_open_HDF5_Chombo(handle, out, opts);
	}

	H5Gclose(gid);
	/*H5Fclose(handle);*/

	return out;
	

}
FR_File *FR_open_HDF5_Pmesh(hid_t handle, FR_File* out, options_t opts){
  int i;
  int status;
  char temp_varnames[FR_MAXVARS*FR_VAR_STRING_SIZE];
  char temp_propnames[FR_MAXVARS*FR_PART_PROP_STRING_SIZE];

  hid_t dataspace, memspace, dataset, name_dataset;
  hsize_t maximum_dims[10];
  hsize_t dataspace_dims[10];
  hsize_t dimens_1d, maxdimens_1d;
  hid_t string_type;
  hid_t int_list_type;
  int_list_t *int_list;

  hid_t num_particles;
  hid_t partTypeId;
  int numElements;
  char *name;
  hid_t type, native_type;
  int numRealPropsFound, numIntPropsFound;

  herr_t (*old_func)(void*);
  void *old_client_data;
  
    int num_vars;
  

  /* integer scalars */
   
  /* grab the data with the name 'integer scalars' from the file */
  dataset = H5Dopen(handle, "integer scalars"); 
  if (dataset < 0) {
    printf("Error opening data file for \"integer scalars\"\n");
    return NULL;
  }

  dataspace = H5Dget_space(dataset);

  /* read the extent of 'dataspace' into 'dimens_1d' */
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);

  /* malloc a pointer to a list of int_list_t's */
  int_list = (int_list_t *) malloc(dimens_1d * sizeof(int_list_t)); 

  /* create an empty vessel sized to hold one int_list_t's worth of data */
  int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
  
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, MAX_STRING_LENGTH);

  /* subdivide the empty vessel into its component sections (name and value) */
  H5Tinsert(int_list_type, 
	    "name", 
	    HOFFSET(int_list_t, name),
	    string_type);

  H5Tinsert(int_list_type, 
	    "value", 
	    HOFFSET(int_list_t, value),
	    H5T_NATIVE_INT);

  /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
  memspace = H5Screate_simple(1, &dimens_1d, NULL);

  status = H5Dread(dataset, int_list_type, memspace, dataspace, H5P_DEFAULT, int_list);
  if (status < 0) {
    printf("Error reading int scalars from data file\n");
    return NULL;
  }

  /* compare this value's 'name' field to the word we're looking for
     using our 'specialcmp' function (defined above) */  
  for (i = 0; i < dimens_1d; i++) {
    if (specialcmp(int_list[i].name, "globalnumblocks") > 0) out->nblocks = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nxb") > 0) out->ncells_vec[0] = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nyb") > 0) out->ncells_vec[1] = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nzb") > 0) out->ncells_vec[2] = int_list[i].value;
  }

  H5Tclose(int_list_type);
  free(int_list);

  H5Sclose(dataspace);
  H5Dclose(dataset);
  /* done with integer scalars */

  /* use ncells_vec to get dim */
  if (out->ncells_vec[2] > 1)
    out->dim = 3;
  else if (out->ncells_vec[1] > 1)
    out->dim = 2;
  else
    out->dim = 1;

  out->ncells = out->ncells_vec[0] *
                out->ncells_vec[1] *
                out->ncells_vec[2] ;

  out->nodetype = (int *) malloc(sizeof(int)*out->nblocks);  
  if(FR_slurp_HDF5(handle, H5T_NATIVE_INT, "node type",out->dim, out->nodetype))
    return NULL;

  out->lref = (int *) malloc(sizeof(int)*out->nblocks);  
  if(FR_slurp_HDF5(handle, H5T_NATIVE_INT, "refine level",out->dim, out->lref))
    return NULL;
  
  out->coord = (double *) malloc(sizeof(double)*out->dim*out->nblocks);
  if(FR_slurp_HDF5(handle, H5T_NATIVE_DOUBLE, "coordinates", out->dim, out->coord))
    return NULL;

  out->size = (double *) malloc(sizeof(double)*out->dim*out->nblocks);
  if(FR_slurp_HDF5(handle, H5T_NATIVE_DOUBLE, "block size", out->dim, out->size))
    return NULL;

  out->bbox = (double *) malloc(sizeof(double)*out->dim*out->nblocks*2);
  if(FR_slurp_HDF5(handle, H5T_NATIVE_DOUBLE, "bounding box", out->dim, out->bbox))
    return NULL;


  /* number of variables and variable names */
  for(i=0; i<(FR_MAXVARS*FR_VAR_STRING_SIZE); i++)
    temp_varnames[i]='\0';

  if(FR_slurp_HDF5(handle, (hid_t) 0, "unknown names", out->dim, temp_varnames))
    return NULL;
  

  if(opts.gridVarSelfDiscovery){
    
    FR_get_data_names_HDF5(handle, &(out->varnames), 
                           &(out->vartypes), &(out->nvar), out->ncells_vec);
  }
  else{ /* this is the way SFOCU originally got grid quantities. */
    out->nvar = 0;
    for(i=0; i<FR_MAXVARS; i++){
      if(temp_varnames[i*FR_VAR_STRING_SIZE]=='\0')
        break;
      (out->nvar)++;
      out->vartypes[i] = UNK;
      strncpy(out->varnames[i], temp_varnames+FR_VAR_STRING_SIZE*i, 
              FR_VAR_STRING_SIZE);
      out->varnames[i][FR_VAR_STRING_SIZE] = '\0';
    }
  }
  FR_GetNumParticles_HDF5(&handle, &num_particles);

  out->totalparticles = num_particles;
  out->numRealPartProps = 0;
  out->numIntPartProps = 0;

  if (out->totalparticles > 0) {
    numRealPropsFound = 0;
    numIntPropsFound = 0;
    /* new file format without compound datatype, and no integer properties */

    /* number of real properties and names */
    for(i=0; i<(FR_MAXVARS*FR_PART_PROP_STRING_SIZE); i++)
      temp_propnames[i]='\0';

    if(FR_slurp_HDF5(handle, (hid_t) 0, "particle names", out->dim, temp_propnames))
      return NULL;

    out->numRealPartProps = 0;
    for(i=0; i<FR_MAXVARS; i++){
      if(temp_propnames[i*FR_PART_PROP_STRING_SIZE]=='\0')
	break;
      (out->numRealPartProps)++;
      strncpy(out->realPartPropNames[i], temp_propnames+FR_PART_PROP_STRING_SIZE*i, 
	      FR_PART_PROP_STRING_SIZE);
      out->realPartPropNames[i][FR_PART_PROP_STRING_SIZE] = '\0';
    }
  }
  return out;  
}

FR_Block *FR_GetBlock_HDF5(FR_File *file, char *var, int block_no){
  FR_Block *out;
  int rank;
  hsize_t dimens_4d[4];
  hssize_t start_4d[5];
  hsize_t stride_4d[5], count_4d[5];
  hid_t ierr;
  herr_t status;
  hid_t dataspace, memspace, dataset;

  int i, varType;
  int xCells, yCells, zCells, numCells;
  
  /*first figure out what kind of block we're looking at*/
  varType = 0;
  for(i=0; i<file->nvar; ++i){
    if(!(strncmp(file->varnames[i], var, FR_VAR_STRING_SIZE))){
      varType = file->vartypes[i];
      break;
    }
  }
  if (!varType) {
    /* try once more with more relaxed comparison... */
    for(i=0; i<file->nvar; ++i){
      if(!(namencmp(file->varnames[i], var, FR_VAR_STRING_SIZE))){
	varType = file->vartypes[i];
	break;
      }
    }
  }

	xCells = file->ncells_vec[0];
	yCells = file->ncells_vec[1];
	zCells = file->ncells_vec[2];
	
	switch(varType){
	case UNK:
		/*do nothing*/  
		break;
	case FACEX:
		xCells += 1;
		break;
	case FACEY:
		yCells += 1;
		break;
	case FACEZ:
		zCells += 1;
		break;
	case SCRATCH:
		xCells += 1;
		yCells += (file->dim >=2)?1:0;
		zCells += (file->dim >=3)?1:0;
		break;
	default:
	  printf("in default case for var \"%s\", variable type unknown: HELP!\n",var);
		/*defaults to unk behavior*/
		break;
		
	}
	numCells = xCells*yCells*zCells;

  out = (FR_Block *) malloc(sizeof(FR_Block));
  out->data = (double *) malloc(sizeof(double)*numCells);

  out->size[0] = xCells;
  out->size[1] = yCells;
  out->size[2] = zCells;

  /* Now do the hdf5 dance */
  
  /* define the dataspace */
  rank = 4;
  dimens_4d[0] = file->nblocks;
  dimens_4d[1] = zCells; /* reverse order */
  dimens_4d[2] = yCells;
  dimens_4d[3] = xCells;

  start_4d[0]  = (hssize_t) block_no;
  start_4d[1]  = 0;
  start_4d[2]  = 0;
  start_4d[3]  = 0;
    
  stride_4d[0] = 1;
  stride_4d[1] = 1;
  stride_4d[2] = 1;
  stride_4d[3] = 1;

  count_4d[0]  = 1; /*Just take one block */
  count_4d[1]  = zCells;
  count_4d[2]  = yCells;
  count_4d[3]  = xCells;

  dataspace = H5Screate_simple(rank, dimens_4d, NULL);
  ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			     start_4d, stride_4d, count_4d, NULL);

  /* define the memspace - just a small modification of the above */
  dimens_4d[0] = 1; /* Again, just one block */
  start_4d[0] = 0;
  memspace = H5Screate_simple(rank, dimens_4d, NULL);
  ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			     start_4d, stride_4d, count_4d, NULL);

  /* Now read */
  dataset = H5Dopen(*((hid_t *)file->handle), var);
  if (dataset < 0) {
    /* Maybe try again, with modified variable name: either more or fewer trailing space characters. */
    char tmpname[6];
    char* p;
    int l;
    strncpy(tmpname, var, 6);
    tmpname[5] = '\0';
    l = strlen(var);
    if (l > 0 && l < 5) {
      if (tmpname[l-1] == ' ') {
	/* eliminate trailing spaces */
	for (l--; (l>0) && ( tmpname[l]==' ');l--) {
	  tmpname[l] = '\0';
	}
      } else {
	/* append spaces */
	for (; (l<4) && (tmpname[l]=='\0');l++) {
	  tmpname[l] = ' ';
	}
      }
    }
    if (strncmp(tmpname,var,6) != 0)
      dataset = H5Dopen(*((hid_t *)file->handle), tmpname);      
  }
  if (dataset < 0) {
    printf("SFOCU has encountered a variable that it cannot read: H5Dopen failed.\nVariable is \"%s\", file \"%s\".\nFAILURE\n",
	   var, file->filename);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    exit(1);
  }
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
		    H5P_DEFAULT, (void *) out->data);
  if (status < 0) {
    printf("SFOCU has encountered a variable that it cannot read: H5Dread failed.\nVariable name is \"%s\".\nFAILURE\n", var);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    exit(1);
  }
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  return out;
}

int FR_get_data_names_HDF5(hid_t filename, char names[FR_MAXVARS][FR_VAR_STRING_SIZE+1], int *var_type, int *num_vars, int *ncells_vec){

  hsize_t num_objs, num_dims, num_obj, obj_type, name_size;
  hid_t dataset, dataspace;

  hsize_t dimens_4d[4], maxdimens_4d[4];
  
  char obj_name[MAX_STRING_LENGTH];

  int i;


  H5Gget_num_objs(filename, &num_objs);

  *num_vars = 0;
  for( i = 0; i < num_objs; ++i){

    /*have to iterate over objects, look at dimensions, and store the ones we want*/
    if(H5Gget_objtype_by_idx(filename, i) != H5G_DATASET)
      continue; /*only interested in datasets*/
    
    H5Gget_objname_by_idx(filename, (hsize_t)i,obj_name, 80);
    
    dataset = H5Dopen(filename, obj_name);
    dataspace = H5Dget_space(dataset);
    if(dataspace < 0)
      printf("dataspace is wrong.\n");
    
    if(H5Sget_simple_extent_ndims(dataspace) != 4){
      H5Sclose(dataspace);
      H5Dclose(dataset);
      continue;
    }
    
    strncpy(names[*num_vars], obj_name, FR_VAR_STRING_SIZE + 1);
        
    /*determine the type of variable we are looking at (face vs unk?)*/
    /*A bit ugly for now, but should work. */
    
     H5Sget_simple_extent_dims(dataspace, dimens_4d, maxdimens_4d);
     
     if(dimens_4d[XDIM] == ncells_vec[0] && 
        dimens_4d[YDIM] == ncells_vec[1] &&
        dimens_4d[ZDIM] == ncells_vec[2])
       var_type[*num_vars] = UNK;
     else if(dimens_4d[XDIM] == ncells_vec[0]+1 && 
        dimens_4d[YDIM] == ncells_vec[1] &&
        dimens_4d[ZDIM] == ncells_vec[2])
       var_type[*num_vars] = FACEX;
     else if(dimens_4d[XDIM] == ncells_vec[0] && 
        dimens_4d[YDIM] == ncells_vec[1]+1 &&
        dimens_4d[ZDIM] == ncells_vec[2])
       var_type[*num_vars] = FACEY;
     else if(dimens_4d[XDIM] == ncells_vec[0] && 
        dimens_4d[YDIM] == ncells_vec[1] &&
        dimens_4d[ZDIM] == ncells_vec[2]+1)
       var_type[*num_vars] = FACEZ;
     else if(dimens_4d[XDIM] == ncells_vec[0]+1 && 
        dimens_4d[YDIM] == ncells_vec[1]+1 &&
        dimens_4d[ZDIM] == ncells_vec[2]+1)
       var_type[*num_vars] = SCRATCH;
     else{
       printf("SFOCU has encountered an variable type it doesn't know how to handle.\nVariable name is \"%s\".\nFAILURE\n", names[*num_vars]);
       exit(1);
     }
              
    H5Sclose(dataspace);
       H5Dclose(dataset);
      (*num_vars)++;
  }
 
  printf("%d\n", *num_vars);
  

  return 0;
}
