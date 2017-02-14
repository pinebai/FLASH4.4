/*get data that lets us figure out how to decompose the domain*/
int FR_getFlashData(char *filename, FR_FileData *fileData){
  
  int i;
  int status;
  hid_t handle;
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

  fileData->handle = malloc(sizeof(hid_t));
  fileData->format = FR_HDF5;
  strcpy(fileData->filename, filename);
  handle =  H5Fopen(fileData->filename, H5F_ACC_READONLY, H5P_DEFAULT);
  *((hid_t*)fileData->handle) = handle;
  
  if (handle < 0)
    return NULL;
  
  
  /* grab the data with the name 'integer scalars' from the file */
  dataset = H5Dopen(handle, "integer scalars"); 
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
    return 1;
  }
    
  /* compare this value's 'name' field to the word we're looking for
     using our 'specialcmp' function (defined above) */  
  for (i = 0; i < dimens_1d; i++) {
    if (specialcmp(int_list[i].name, "globalnumblocks") > 0) flashData->nblocks = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nxb") > 0) flashData->ncells_vec[0] = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nyb") > 0) flashData->ncells_vec[1] = int_list[i].value;
    else if (specialcmp(int_list[i].name, "nzb") > 0) flashData->ncells_vec[2] = int_list[i].value;
  }

  H5Tclose(int_list_type);
  free(int_list);

  H5Sclose(dataspace);
  H5Dclose(dataset);
  /* done with integer scalars */

  /* use ncells_vec to get dim */
  if (flashData->ncells_vec[2] > 1)
    flashData->dim = 3;
  else if (flashData->ncells_vec[1] > 1)
    flashData->dim = 2;
  else
    flashData->dim = 1;

  flashData->ncells = flashData->ncells_vec[0] *
                flashData->ncells_vec[1] *
                flashData->ncells_vec[2] ;


  /* number of variables and variable names */
  for(i=0; i<(FR_MAXVARS*FR_VAR_STRING_SIZE); i++)
    temp_varnames[i]='\0';

  if(FR_slurp_HDF5(handle, (hid_t) 0, "unknown names", temp_varnames))
    return NULL;
  
  flashData->nvar = 0;
  for(i=0; i<FR_MAXVARS; i++){
    if(temp_varnames[i*FR_VAR_STRING_SIZE]=='\0')
      break;
    (flashData->nvar)++;
    strncpy(flashData->varnames[i], temp_varnames+FR_VAR_STRING_SIZE*i, 
	    FR_VAR_STRING_SIZE);
    flashData->varnames[i][FR_VAR_STRING_SIZE] = '\0';
  }

  FR_GetNumParticles_HDF5(&handle, &num_particles);

  flashData->totalparticles = num_particles;
  flashData->numRealPartProps = 0;
  flashData->numIntPartProps = 0;

  if (flashData->totalparticles > 0) {
    numRealPropsFound = 0;
    numIntPropsFound = 0;
    /* new file format without compound datatype, and no integer properties */

    /* number of real properties and names */
    for(i=0; i<(FR_MAXVARS*FR_PART_PROP_STRING_SIZE); i++)
      temp_propnames[i]='\0';

    if(FR_slurp_HDF5(handle, (hid_t) 0, "particle names", temp_propnames))
      return NULL;

    flashData->numRealPartProps = 0;
    for(i=0; i<FR_MAXVARS; i++){
      if(temp_propnames[i*FR_PART_PROP_STRING_SIZE]=='\0')
	break;
      (flashData->numRealPartProps)++;
      strncpy(flashData->realPartPropNames[i], temp_propnames+FR_PART_PROP_STRING_SIZE*i, 
	      FR_PART_PROP_STRING_SIZE);
      flashData->realPartPropNames[i][FR_PART_PROP_STRING_SIZE] = '\0';
    }
  }

  return 0;
}
