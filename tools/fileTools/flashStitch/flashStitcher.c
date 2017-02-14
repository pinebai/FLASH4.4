#include "flashStitcher.h"
#include <string.h>
#include "hdf5.h"

/*read in a group of output files and output as one hdf5 file*/

/*USAGE: flashStitcher [options] finalFilename 
  final file name is assumed to match format of the input file names for now*/

int flashStitcher(char *filename, int maxSplitNum, int splitZeroPad,
                  options_t *opts){

  int i; /*counters.  We're gonna need these*/

  /*hdf5 handles, etc*/
  herr_t status;
  hsize_t dataspace_dims[10];
  hid_t outFile, inFile;

  int currentSplitNum = 0;
  
  /*filename formatting.  let's make finding the next file easy.*/
  char splitFilenameFormat[STRINGSIZE];
  char splitFilename [STRINGSIZE];
  char splitNumFmt[20];
  char filenameTokens[STRINGSIZE];
  char *basenm; 
  char *storageFormat;
  char *fileType;
  char *fileNum;


  /*hdf5 handles*/
  hid_t dataspaceIn, dataspaceOut, datasetIn, datasetOut;
  hid_t memspaceIn, memspaceOut;
  hsize_t dimens_1d, maxdimens_1d;

  hid_t real_list_type, int_list_type, str_list_type, log_list_type;
  hid_t string_type;
  
  int_list_t *int_list;
  real_list_t *real_list;
  str_list_t *str_list;
  log_list_t *log_list;
  sim_info_t simInfo;
  
  int numUnkNames;
  char *unkNames;

  int globalNumBlocks;
  int currentOffset = 0;
  int nextOffset = 0;
  int firstCall = 1;

  /*figure out the split filenames from the passed in name.
    overridable later...*/
  /* printf("starting to build filename format %s\n", filename);
  strncpy(filenameTokens, filename, (size_t)STRINGSIZE);

  printf("cpy good.\n");
  

  basenm = strtok(filenameTokens, "_");
  printf("basenm = %s\n", basenm);

  storageFormat = strtok(NULL, "_");
  printf("storageFormat = %s\n", storageFormat);

  fileType = strtok(NULL, "_");
  printf("fileType = %s %d\n", fileType, strncmp(fileType, "chk", 3));

  if(0 != strncmp(fileType, "chk", 3)){
    printf("why are we here?\n");
    fflush(stdout);
    strtok(NULL, "_");
    strncat(fileType,"_cnt", STRINGSIZE);
  }
  

  fileNum = strtok(NULL, "_");
  printf("fileNum = %s\n", fileNum);
  */

  printf("constructing format\n");
  sprintf(splitNumFmt, "s%%0%dd", splitZeroPad);
    printf("splitNumFmt = %s\n", splitNumFmt);  
    sprintf(splitFilenameFormat, "%s_%s_%s_%s_%s", opts->basenm,splitNumFmt,
            opts->format, opts->type, opts->filenum); 
    /* storageFormat, fileType, fileNum);*/
      printf("splitFilenameFormat = %s\n", splitFilenameFormat);
  
 
  printf("%s\n%s\n%s\n%s\n%s\n%s\n", basenm, splitNumFmt, storageFormat, fileType, fileNum, splitFilenameFormat);
 /*open out destination file*/

  sprintf(splitFilename, splitFilenameFormat, 3);
  printf("simulated filename = %s\n", splitFilename);
  sprintf(splitFilename, splitFilenameFormat, 4);
  printf("simulated filename = %s\n", splitFilename);

  printf("opening %s\n", filename);
  outFile = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT ,H5P_DEFAULT);
  printf("outputFile opened!\n");


  /*iterate over split files and copy them*/ 
  while(currentSplitNum < maxSplitNum){
    
    sprintf(splitFilename, splitFilenameFormat, currentSplitNum);
    inFile = H5Fopen(splitFilename,H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if(currentSplitNum == 0){ /*get metadata and copy it*/
      /* Yes, this is ugly.  I havne't thought of a better way to do this
         yet, though.*/
      
      /*stuff that is kind of universal*/
      /*out general string type*/
      string_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(string_type, MAX_STRING_LENGTH);
      
      /*integer scalars******************************************************/
      datasetIn = H5Dopen(inFile, "integer scalars");
      
      /*make sure to grab sizes that we need from here now.*/
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);
      
      int_list = (int_list_t*)malloc(dimens_1d *sizeof(int_list_t));
      
      /*create datatype ++++*/
      int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
      
      H5Tinsert(int_list_type, "name", HOFFSET(int_list_t, name), string_type);
      H5Tinsert(int_list_type, "value", HOFFSET(int_list_t, value), H5T_NATIVE_INT);
      
      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);
      
      status = H5Dread(datasetIn, int_list_type, memspaceIn, dataspaceIn, 
                       H5P_DEFAULT, int_list);
      if(status <0){
        /*we should fail gracefully*/
      }
      
      /*write this back out*/
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "integer scalars", int_list_type,
                             dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, int_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        int_list);
      if(status <0){
        /*more failing gracefuly*/
      }
      
      /*we need a couple of values from here*/
      for(i = 0; i < dimens_1d; ++i){
        if(0 == strncmp(int_list[i].name, "globalnumblocks",15))
          globalNumBlocks = int_list[i].value;
      }
 
      /*close things up for next iteration.*/
      free(int_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);

      /*integer runtime parameters*****************************************/
      
      datasetIn = H5Dopen(inFile, "integer runtime parameters");
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);
      
      int_list = (int_list_t *) malloc(dimens_1d*sizeof(int_list_t));
      
      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);
     
      status = H5Dread(datasetIn, int_list_type, memspaceIn, dataspaceIn,
                       H5P_DEFAULT, int_list);

      if(status <0){
        /*hey, fail gracefully agian!*/
      }
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "integer runtime parameters", 
                             int_list_type, dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, int_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        int_list);
      
      if(status <0){
        /*Guess, what? We FAILED again*/
      }

      /*close things up for next iteration.*/
      free(int_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);

      /*clean up types*/
      H5Tclose(int_list_type);



      /*real scalars******************************************************/
      datasetIn = H5Dopen(inFile, "real scalars");

      /*make sure to grab sizes that we need from here now.*/
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);

      real_list = (real_list_t*) malloc(dimens_1d *sizeof(real_list_t));
     
      /*create datatype ++++*/
      real_list_type = H5Tcreate(H5T_COMPOUND, sizeof(real_list_t));

      H5Tinsert(real_list_type, "name", HOFFSET(real_list_t, name), string_type);
      H5Tinsert(real_list_type, "value", HOFFSET(real_list_t, value), H5T_NATIVE_DOUBLE);

      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);

      status = H5Dread(datasetIn, real_list_type, memspaceIn, dataspaceIn, 
                       H5P_DEFAULT, real_list);
      if(status <0){
        /*we should fail gracefully*/
      }

      /*write this back out*/
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "real scalars", real_list_type,
                             dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, real_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        real_list);
      if(status <0){
        /*more failing gracefuly*/
      }
      
      /*close things up for next iteration.*/
      free(real_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);
      
      /*real runtime parameters*****************************************/
      
      datasetIn = H5Dopen(inFile, "real runtime parameters");
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);
      
      real_list = (real_list_t *) malloc(dimens_1d*sizeof(real_list_t));
      
      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);
     
      status = H5Dread(datasetIn, real_list_type, memspaceIn, dataspaceIn,
                       H5P_DEFAULT, real_list);

      if(status <0){
        /*hey, fail gracefully agian!*/
      }
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "real runtime parameters", 
                             real_list_type, dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, real_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        real_list);
      
      if(status <0){
        /*Guess, what? We FAILED again*/
      }

      /*close things up for next iteration.*/
      free(real_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);

      /*clean up types*/
      H5Tclose(real_list_type);

      /*string scalars******************************************************/
      datasetIn = H5Dopen(inFile, "string scalars");

      /*make sure to grab sizes that we need from here now.*/
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);

      str_list = (str_list_t*) malloc(dimens_1d *sizeof(str_list_t));
     
      /*create datatype ++++*/
      str_list_type = H5Tcreate(H5T_COMPOUND, sizeof(str_list_t));

      H5Tinsert(str_list_type, "name", HOFFSET(str_list_t, name), 
                string_type);
      H5Tinsert(str_list_type, "value", HOFFSET(str_list_t, value), 
                string_type);

      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);

      status = H5Dread(datasetIn, str_list_type, memspaceIn, dataspaceIn, 
                       H5P_DEFAULT, str_list);
      if(status <0){
        /*we should fail gracefully*/
      }

      /*write this back out*/
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "string scalars", str_list_type,
                             dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, str_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        str_list);
      if(status <0){
        /*more failing gracefuly*/
      }
        
      /*close things up for next iteration.*/
      free(str_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);
      
      /*string runtime parameters*****************************************/
      
      datasetIn = H5Dopen(inFile, "string runtime parameters");
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);
      
      str_list = (str_list_t *) malloc(dimens_1d*sizeof(str_list_t));
      
      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);
     
      status = H5Dread(datasetIn, str_list_type, memspaceIn, dataspaceIn,
                       H5P_DEFAULT, str_list);

      if(status <0){
        /*hey, fail gracefully agian!*/
      }
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "string runtime parameters", 
                             str_list_type, dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, str_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        str_list);
      
      if(status <0){
        /*Guess, what? We FAILED again*/
      }

      /*close things up for next iteration.*/
      free(str_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);

      /*clean up types*/
      H5Tclose(str_list_type);



      /*logical scalars******************************************************/
      datasetIn = H5Dopen(inFile, "logical scalars");

      /*make sure to grab sizes that we need from here now.*/
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);

      log_list = (log_list_t*) malloc(dimens_1d *sizeof(log_list_t));
     
      /*create datatype ++++*/
      log_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));

      H5Tinsert(log_list_type, "name", HOFFSET(log_list_t, name), string_type);
      H5Tinsert(log_list_type, "value", HOFFSET(log_list_t, value), 
                H5T_NATIVE_INT);

      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);

      status = H5Dread(datasetIn, log_list_type, memspaceIn, dataspaceIn, 
                       H5P_DEFAULT, log_list);
      if(status <0){
        /*we should fail gracefully*/
      }

      /*write this back out*/
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "logical scalars", log_list_type,
                             dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, log_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        log_list);
      if(status <0){
        /*more failing gracefuly*/
      }
        
      /*close things up for next iteration.*/
      free(log_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);
      
      /*logical runtime parameters*****************************************/
      
      datasetIn = H5Dopen(inFile, "logical runtime parameters");
      dataspaceIn = H5Dget_space(datasetIn);
      H5Sget_simple_extent_dims(dataspaceIn, &dimens_1d, &maxdimens_1d);
      
      log_list = (log_list_t *) malloc(dimens_1d*sizeof(log_list_t));
      
      memspaceIn = H5Screate_simple(1, &dimens_1d, NULL);
     
      status = H5Dread(datasetIn, log_list_type, memspaceIn, dataspaceIn,
                       H5P_DEFAULT, log_list);

      if(status <0){
        /*hey, fail gracefully agian!*/
      }
      dataspaceOut = H5Screate_simple(1, &dimens_1d, NULL);
      datasetOut = H5Dcreate(outFile, "logical runtime parameters", 
                             log_list_type, dataspaceOut, H5P_DEFAULT);
      status = H5Dwrite(datasetOut, log_list_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        log_list);
      
      if(status <0){
        /*Guess, what? We FAILED again*/
      }

      /*close things up for next iteration.*/
      free(log_list);
      
      H5Dclose(datasetIn);
      H5Sclose(dataspaceIn);
      H5Sclose(memspaceIn);
      H5Dclose(datasetOut);
      H5Sclose(dataspaceOut);

      /*clean up types*/
      H5Tclose(log_list_type);



      /*Now copy header metadata.********************************************/
      /*we also get the unk labels from here*/
      /*see if we have particles, and if so, get their names, too*/

      /*DEV: TODO:will do particles later*/
      
      readFlashInfo(&inFile, &simInfo, &unkNames, &numUnkNames);
      printSimInfo(&simInfo);
      printf("RETURNED SUCCESSFULLY\n");
      writeFlashInfo(&outFile, &simInfo, unkNames, &numUnkNames);
      

    } /*end metadata copy*/
    

    /***********************************************************************/
     /*all files participate from here on out*/

    /*copy bflags*/
    copy_bflags(inFile, outFile, currentOffset, globalNumBlocks, &nextOffset);
    printf("CURRENT OFFSET: %d NEXT OFFSET: %d\n", currentOffset, nextOffset);
    /*copy gid*/
    copy_gid(inFile, outFile, currentOffset, globalNumBlocks);
    /*copy block size*/
    copy_blocksize(inFile, outFile, currentOffset, globalNumBlocks);
    /*copy bound box*/
    copy_bndbox(inFile, outFile, currentOffset, globalNumBlocks);
    /*copy coordinates*/
    copy_coordinates(inFile, outFile, currentOffset, globalNumBlocks);
    /*copy node type*/
    copy_1d_dataset(inFile, outFile, currentOffset, globalNumBlocks,
                    H5T_NATIVE_INT, "node type", firstCall);
    /*copy processor number*/
    copy_1d_dataset(inFile, outFile, currentOffset, globalNumBlocks,
                    H5T_NATIVE_INT, "processor number", firstCall);
    /*copy refine level*/
    copy_1d_dataset(inFile, outFile, currentOffset, globalNumBlocks,
                    H5T_NATIVE_INT, "refine level",firstCall);
    /*copy which child*/
    copy_1d_dataset(inFile, outFile, currentOffset, globalNumBlocks,
                    H5T_NATIVE_INT, "which child",firstCall);

    printf("STARTING ON UNKVARS\n");
    /*copy all unknown variables*/
    copy_unkvars(inFile, outFile, currentOffset, globalNumBlocks,
                 unkNames, numUnkNames, firstCall);
    


    
    firstCall = 0;
    currentOffset = nextOffset;
    H5Fclose(inFile);
    /*break;*/
    currentSplitNum++;
    
  }
  

  /* deallocate unkNames*/
  free(unkNames);
  H5Fclose(outFile);
  
  return 0;
  
}


int readFlashInfo(hid_t *file,sim_info_t *siminfo, char **unkNames, int *numUnkNames){
  
  
  hid_t dataset, dataspace, memspace;
  herr_t status;
  hsize_t dimens_1d, dimens_2d[2], maxdimens_1d, maxdimens_2d[2];
  
  int string_size;
  
  hid_t string_type, sp_type, string_type_setup, string_type_max_str_len, si_type;
  
  hid_t attribute, attribute_space;
  int rank;
  hid_t group_id;
  int i = 0;
  
  /*build siminfo datatype*/
  string_type_max_str_len = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_max_str_len, MAX_STRING_LENGTH);
  
  string_type_setup = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_setup, 400);
  
  rank = 1;
  
  si_type = H5Tcreate(H5T_COMPOUND, sizeof(sim_info_t));
    
  H5Tinsert(si_type, 
            "file format version", 
            offsetof(sim_info_t, file_format_version),
            H5T_NATIVE_INT);
  
  H5Tinsert(si_type, 
            "setup call", 
            HOFFSET(sim_info_t, setup_call),
            string_type_setup);
  
  H5Tinsert(si_type, 
            "file creation time", 
            HOFFSET(sim_info_t, file_creation_time),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "flash version", 
            HOFFSET(sim_info_t, flash_version),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build date", 
            HOFFSET(sim_info_t, build_date),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build dir", 
            HOFFSET(sim_info_t, build_dir),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build machine", 
            HOFFSET(sim_info_t, build_machine),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "cflags", 
            HOFFSET(sim_info_t, cflags),
            string_type_setup);
  
  
  H5Tinsert(si_type, 
            "fflags", 
            HOFFSET(sim_info_t, fflags),
            string_type_setup);
  
  
  H5Tinsert(si_type, 
            "setup time stamp", 
            HOFFSET(sim_info_t,setup_time_stamp),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build time stamp", 
            HOFFSET(sim_info_t, build_time_stamp),
            string_type_max_str_len);
  
  /*read in dataset*/
  dataset = H5Dopen(*file, "sim info");
  dataspace = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
  
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  
  status = H5Dread(dataset, si_type, memspace, dataspace, H5P_DEFAULT, siminfo);
  printf("SimInfo Read!\n");
  fflush(stdout);
  
  H5Tclose(string_type_setup);
  H5Tclose(string_type_max_str_len);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(si_type);
  
  
  /*now output unknown names*/
  rank  = 2;
  dimens_2d[0] = (hsize_t) *numUnkNames;
  dimens_2d[1] = 1;
  
  dataset = H5Dopen(*file, "unknown names");
  dataspace = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
  *numUnkNames = dimens_2d[0];

  
  *unkNames = (char *)malloc((*numUnkNames)*UNK_NAME_LEN*sizeof(char*));

  /*still expecting that string size of 4...*/
  string_size = 4;
  
  string_type = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(string_type, string_size);
  
  memspace = H5Screate_simple(rank, dimens_2d, NULL);

  status = H5Dread(dataset, string_type, memspace, dataspace, H5P_DEFAULT, *unkNames);  
  for(i = 0; i < 10; ++i){
  }
  
  printf("unkNames read!\n");  
  H5Tclose(string_type);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  printf("about to return\n");
  
  return 0;
}

int writeFlashInfo(hid_t *file, sim_info_t *siminfo, char *unkNames, int *numUnkNames){

  hid_t dataset, dataspace, memspace;
  herr_t status;
  hsize_t dimens_1d, dimens_2d[2];

  int string_size;

  hid_t string_type, sp_type, string_type_setup, string_type_max_str_len, si_type;

  hid_t attribute, attribute_space;
  int rank;
  hid_t group_id;
  int i = 0;


  /*build siminfo datatype*/
  string_type_max_str_len = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_max_str_len, MAX_STRING_LENGTH);



  string_type_setup = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_setup, 400);

  rank = 1;
  dimens_1d = 1;
  
  dataspace = H5Screate_simple(rank, &dimens_1d, NULL);
  si_type = H5Tcreate(H5T_COMPOUND, sizeof(sim_info_t));
  
  H5Tinsert(si_type, 
            "file format version", 
            offsetof(sim_info_t, file_format_version),
            H5T_NATIVE_INT);
  
  H5Tinsert(si_type, 
            "setup call", 
            HOFFSET(sim_info_t, setup_call),
            string_type_setup);
  
  H5Tinsert(si_type, 
            "file creation time", 
            HOFFSET(sim_info_t, file_creation_time),
            string_type_max_str_len);
  
  H5Tinsert(si_type,
            "flash version", 
            HOFFSET(sim_info_t, flash_version),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build date", 
            HOFFSET(sim_info_t, build_date),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build dir", 
            HOFFSET(sim_info_t, build_dir),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build machine", 
            HOFFSET(sim_info_t, build_machine),
            string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "cflags", 
            HOFFSET(sim_info_t, cflags),
            string_type_setup);
  
  
  H5Tinsert(si_type, 
            "fflags", 
            HOFFSET(sim_info_t, fflags),
            string_type_setup);
  
  
  H5Tinsert(si_type,
            "setup time stamp", 
            HOFFSET(sim_info_t,setup_time_stamp),
                        string_type_max_str_len);
  
  H5Tinsert(si_type, 
            "build time stamp", 
            HOFFSET(sim_info_t, build_time_stamp),
            string_type_max_str_len);
  
  
  
  dataset = H5Dcreate(*file, "sim info", si_type,
                      dataspace, H5P_DEFAULT);
  
  status = H5Dwrite(dataset, si_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, siminfo);
  
  
  H5Tclose(string_type_setup);
  H5Tclose(string_type_max_str_len);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(si_type);
  
  
  printf("starting unknames\n");
  /*now output unknown names*/
  rank  = 2;
  dimens_2d[0] = (hsize_t) *numUnkNames;
  dimens_2d[1] = 1;
  
  
  /*still expecting that string size of 4...*/
  string_size = 4;
  
  string_type = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(string_type, string_size);
  
  dataspace = H5Screate_simple(rank, dimens_2d, NULL);
  dataset   = H5Dcreate(*file, "unknown names", 
                        string_type, dataspace, H5P_DEFAULT);
  printf("copying UnkNames\n");

  
  status = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, unkNames);
  printf("unkNmaes copied\n");
  H5Tclose(string_type);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  return 0;
}

/* utility print function*/
void printSimInfo(sim_info_t *si){

  printf("file_format_version: %d\n", si->file_format_version);
  printf("setup call:\n%s\n", si->setup_call);
  printf("file_creation_time: %s\n", si->file_creation_time);
  printf("flash_version: %s\n", si->flash_version);
  printf("build_date: %s\n", si->build_date);
  printf("build_dir: %s\n", si->build_dir);
  printf("build_machine: %s\n", si->build_machine);
  printf("cflags: %s\n", si->cflags);
  printf("fflags: %s\n", si->fflags);
  printf("setup_time_stamp: %s\n", si->setup_time_stamp);
  printf("build_time_stamp: %s\n", si->build_time_stamp);
  return;
}

/*copy bflags*/

void copy_bflags(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks, int *nextOffset){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hsize_t start_2d[2], count_2d[2], stride_2d[2];
  hsize_t fulldimens_2d[2];

  int *bflags;
 
  static int firstCall = 1;


  rank = 2;
  /* open the dataset */
  dataset = H5Dopen(infile, "bflags");

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);

   
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  /*figure out how far we go in this file*/
  *nextOffset += dimens_2d[0];
  printf("nextOffset = %d\n", *nextOffset);
  
  bflags = (int*)malloc(sizeof(int)*dimens_2d[0]*dimens_2d[1]);
  
  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                H5P_DEFAULT, bflags);

  if (status < 0){
    printf("Error: Unable to read dataset for bflags\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("bflags read in!\n");
  fflush(stdout);
  rank = 2;

  start_2d[0] = (hsize_t)currentOffset;
  start_2d[1] =0;
  
  stride_2d[0] = 1;
  stride_2d[1] = 1;  

  count_2d[0] = dimens_2d[0];
  count_2d[1] = dimens_2d[1];
  
  fulldimens_2d[0] = (hsize_t)globalNumBlocks;
  fulldimens_2d[1] = dimens_2d[1];
  dataspace = H5Screate_simple(rank, fulldimens_2d, NULL);
  
  if(firstCall){
    firstCall = 0;
    printf("creating...\n");
    dataset = H5Dcreate(outfile, "bflags", H5T_NATIVE_INT, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, "bflags");
  }

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                               stride_2d, count_2d, NULL);
  
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  
  status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace,dataspace, 
                    H5P_DEFAULT, bflags);

  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(bflags);

  return;
}
/*copy blocksize*/
 

void copy_blocksize(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hsize_t start_2d[2], count_2d[2], stride_2d[2];
  hsize_t fulldimens_2d[2];

  double *blksize;

  static int firstCall = 1;


  rank = 2;
  /* open the dataset */
  dataset = H5Dopen(infile, "block size");

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);

   
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  /*figure out how far we go in this file*/
   
  blksize = (double*)malloc(sizeof(double)*dimens_2d[0]*dimens_2d[1]);
  
  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                H5P_DEFAULT, blksize);

  if (status < 0){
    printf("Error: Unable to read dataset for blocksize\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("blocksize read in!\n");
  fflush(stdout);
  rank = 2;

  start_2d[0] = (hsize_t)currentOffset;
  start_2d[1] =0;
  
  stride_2d[0] = 1;
  stride_2d[1] = 1;  

  count_2d[0] = dimens_2d[0];
  count_2d[1] = dimens_2d[1];
  
  fulldimens_2d[0] = (hsize_t)globalNumBlocks;
  fulldimens_2d[1] = dimens_2d[1];
  dataspace = H5Screate_simple(rank, fulldimens_2d, NULL);
  
  if(firstCall){
    firstCall = 0;
    dataset = H5Dcreate(outfile, "block size", H5T_NATIVE_DOUBLE, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, "block size");
  }

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                               stride_2d, count_2d, NULL);
  
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,dataspace, 
                    H5P_DEFAULT, blksize);

  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(blksize);

  return;
}

/*copy coordinates*/
void copy_coordinates(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hsize_t start_2d[2], count_2d[2], stride_2d[2];
  hsize_t fulldimens_2d[2];

  double *coords;

  static int firstCall = 1;


  rank = 2;
  /* open the dataset */
  dataset = H5Dopen(infile, "coordinates");

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);

   
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  /*figure out how far we go in this file*/
   
  coords = (double*)malloc(sizeof(double)*dimens_2d[0]*dimens_2d[1]);
  
  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                H5P_DEFAULT, coords);
  if (status < 0){
    printf("Error: Unable to read dataset for coordinates.\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("coordinates read in!\n");
  fflush(stdout);
  rank = 2;

  start_2d[0] = (hsize_t)currentOffset;
  start_2d[1] =0;
  
  stride_2d[0] = 1;
  stride_2d[1] = 1;  

  count_2d[0] = dimens_2d[0];
  count_2d[1] = dimens_2d[1];
  
  fulldimens_2d[0] = (hsize_t)globalNumBlocks;
  fulldimens_2d[1] = dimens_2d[1];
  dataspace = H5Screate_simple(rank, fulldimens_2d, NULL);
  
  if(firstCall){
    firstCall = 0;
    dataset = H5Dcreate(outfile, "coordinates", H5T_NATIVE_DOUBLE, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, "coordinates");
  }

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                               stride_2d, count_2d, NULL);
  
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,dataspace, 
                    H5P_DEFAULT, coords);

  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(coords);

  return;
}

/*copy gid*/
void copy_gid(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hsize_t start_2d[2], count_2d[2], stride_2d[2];
  hsize_t fulldimens_2d[2];

  int *gid;

  static int firstCall = 1;


  rank = 2;
  /* open the dataset */
  dataset = H5Dopen(infile, "gid");

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);

   
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  /*figure out how far we go in this file*/
   
  gid = (int*)malloc(sizeof(int)*dimens_2d[0]*dimens_2d[1]);
  
  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                H5P_DEFAULT, gid);
  if (status < 0){
    printf("Error: Unable to read dataset for gid.\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("gid read in!\n");
  fflush(stdout);
  rank = 2;

  start_2d[0] = (hsize_t)currentOffset;
  start_2d[1] =0;
  
  stride_2d[0] = 1;
  stride_2d[1] = 1;  

  count_2d[0] = dimens_2d[0];
  count_2d[1] = dimens_2d[1];
  
  fulldimens_2d[0] = (hsize_t)globalNumBlocks;
  fulldimens_2d[1] = dimens_2d[1];
  dataspace = H5Screate_simple(rank, fulldimens_2d, NULL);
  
  if(firstCall){
    firstCall = 0;
    dataset = H5Dcreate(outfile, "gid", H5T_NATIVE_INT, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, "gid");
  }

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                               stride_2d, count_2d, NULL);
  
  memspace = H5Screate_simple(rank, dimens_2d, NULL);
  
  status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace,dataspace, 
                    H5P_DEFAULT, gid);

  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(gid);

  return;
}


/*copy bounding box*/
void copy_bndbox(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_3d[3], maxdimens_3d[3];
  hsize_t start_3d[3], count_3d[3], stride_3d[3];
  hsize_t fulldimens_3d[3];

  double *bndbox;

  static int firstCall = 1;

   rank = 3;
  /* open the dataset */
  dataset = H5Dopen(infile, "bounding box");

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_3d, maxdimens_3d);

   
  memspace = H5Screate_simple(rank, dimens_3d, NULL);
  /*figure out how far we go in this file*/
   
  bndbox = (double*)malloc(sizeof(double)*
                           dimens_3d[0]*dimens_3d[1]*dimens_3d[2]);
  
  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                H5P_DEFAULT, bndbox);
  if (status < 0){
    printf("Error: Unable to read dataset for bounding box.\n");
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("bounding box read in!\n");
  fflush(stdout);
  rank = 3;
 printf("currentOffset: %d\n", currentOffset);

  start_3d[0] = (hsize_t)currentOffset;
  start_3d[1] = 0;
  start_3d[2] = 0;
  
  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;  

  printf("current count is: %d\n offset is %d\n", dimens_3d[0], start_3d[0]);
  count_3d[0] = dimens_3d[0];
  count_3d[1] = dimens_3d[1];
  count_3d[2] = dimens_3d[2];


  fulldimens_3d[0] = (hsize_t)globalNumBlocks;
  fulldimens_3d[1] = dimens_3d[1];
  fulldimens_3d[2] = dimens_3d[2];

  dataspace = H5Screate_simple(rank, fulldimens_3d, NULL);
  printf("dataspace created.\n");
  if(firstCall){
    firstCall = 0;
    dataset = H5Dcreate(outfile, "bounding box", H5T_NATIVE_DOUBLE, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, "bounding box");
  }

  printf("dataset open.\n");
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d,
                               stride_3d, count_3d, NULL);
  printf("hyperslab selected]\n");
  memspace = H5Screate_simple(rank, dimens_3d, NULL);
  printf("memspace done.\n");
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,dataspace, 
                    H5P_DEFAULT, bndbox);
  printf("write complete.\n");
  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(bndbox);

  return;
}


/*copy generic nblocks dataspace*/

void copy_1d_dataset(hid_t infile, hid_t outfile, int currentOffset, 
                     int globalNumBlocks, hid_t type, char *name, 
                     int firstCall){
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_1d, maxdimens_1d;
  hsize_t start_1d, count_1d, stride_1d;
  hsize_t fulldimens_1d;

  void *data;


   rank = 1;
  /* open the dataset */
  dataset = H5Dopen(infile, name);

  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);

   
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  /*figure out how far we go in this file*/
  

  if(type == H5T_NATIVE_INT){
    data = malloc(sizeof(int)*(int)dimens_1d);
  } 
  else if(type == H5T_NATIVE_DOUBLE){
    data = malloc(sizeof(double)*(int)dimens_1d);
  }
  
  /* read the data */
  status = H5Dread(dataset, type, memspace, dataspace, 
                H5P_DEFAULT, data);
  if (status < 0){
    printf("Error: Unable to read dataset for %s\n", name);
  }

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  printf("%s read in!\n",name);
  fflush(stdout);
  rank = 1;
 
  start_1d = (hsize_t)currentOffset;
    
  stride_1d = 1;
 
  count_1d = dimens_1d;
 


  fulldimens_1d = (hsize_t)globalNumBlocks;
  
  dataspace = H5Screate_simple(rank, &fulldimens_1d, NULL);
  printf("dataspace created.\n");
  if(firstCall){
    dataset = H5Dcreate(outfile, name, type, 
                        dataspace, H5P_DEFAULT);
  }
  else{
    dataset = H5Dopen(outfile, name);
  }

  printf("dataset open.\n");
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d,
                               &stride_1d, &count_1d, NULL);
  printf("hyperslab selected]\n");
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  printf("memspace done.\n");
  status = H5Dwrite(dataset, type, memspace,dataspace, 
                    H5P_DEFAULT, data);
  printf("write complete.\n");
  H5Sclose(memspace);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  free(data);

  return;
}

/*copy unk*/
/*copy bounding box*/
void copy_unkvars(hid_t infile, hid_t outfile, int currentOffset, 
                  int globalNumBlocks, char* unkNames, int numUnkNames,
                  int firstCall){
  hid_t dataspace, dataset, memspace, attribute_space, attribute;
  herr_t status;

  int rank;
  hsize_t dimens_4d[4], maxdimens_4d[4];
  hsize_t start_4d[4], count_4d[4], stride_4d[4];
  hsize_t fulldimens_4d[4];
  hsize_t dimens_1d;

  double *unk;
  double max, min;
  char name[5];
  int i;
  
  for(i = 0; i < numUnkNames; ++i){
    /*strncpy(name, unkNames[i*4],4);*/
    name[0] = unkNames[i*4];
    name[1] = unkNames[i*4+1];
    name[2] = unkNames[i*4+2];
    name[3] = unkNames[i*4+3];
    name[4] = '\0';
    printf("%s\n", name);
    
    
    /*open for reading*/
    rank = 4;
    
    dataset = H5Dopen(infile, name);
    
    dataspace = H5Dget_space(dataset);
    
    H5Sget_simple_extent_dims(dataspace, dimens_4d, maxdimens_4d);
    
    memspace = H5Screate_simple(rank, dimens_4d, NULL);
    
    
    unk = (double*)malloc(sizeof(double)*
                          dimens_4d[0]*dimens_4d[1]*dimens_4d[2]*dimens_4d[3]);
    
    /*get min and max from first file*/
    if(firstCall){
      attribute = H5Aopen_name(dataset, "minimum");
      H5Aread(attribute, H5T_NATIVE_DOUBLE, &min);
      printf("MIN: %f\n", min);
      H5Aclose(attribute);
      
      attribute = H5Aopen_name(dataset, "maximum");
      H5Aread(attribute, H5T_NATIVE_DOUBLE, &max);
      printf("MAX: %f\n", max);
      H5Aclose(attribute);
    }

    
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                     H5P_DEFAULT, unk);
    if (status < 0){
      printf("Error: Unable to read dataset for %s.\n", name);
    }
    
    H5Sclose(memspace); 
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    printf("%s read in!\n", name);
    fflush(stdout);
    rank = 4;
    printf("currentOffset: %d\n", currentOffset);
    
    start_4d[0] = (hsize_t)currentOffset;
    start_4d[1] = 0;
    start_4d[2] = 0;
    start_4d[3] = 0;
    
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;  
    stride_4d[3] = 1;
    
    
    printf("current count is: %d\n offset is %d\n", dimens_4d[0], start_4d[0]);
    count_4d[0] = dimens_4d[0];
    count_4d[1] = dimens_4d[1];
    count_4d[2] = dimens_4d[2];
    count_4d[3] = dimens_4d[3];
    
    fulldimens_4d[0] = (hsize_t)globalNumBlocks;
    fulldimens_4d[1] = dimens_4d[1];
    fulldimens_4d[2] = dimens_4d[2];
    fulldimens_4d[3] = dimens_4d[3];
    
    dataspace = H5Screate_simple(rank, fulldimens_4d, NULL);
    printf("dataspace created.\n");
    if(firstCall){
      dataset = H5Dcreate(outfile, name, H5T_NATIVE_DOUBLE, 
                          dataspace, H5P_DEFAULT);
      dimens_1d = 1;
      attribute_space = H5Screate_simple(1,&dimens_1d, NULL);
      attribute = H5Acreate(dataset, "minimum", H5T_NATIVE_DOUBLE,
                            attribute_space, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_DOUBLE, &min);
      H5Aclose(attribute);
      attribute = H5Acreate(dataset,"maximum", H5T_NATIVE_DOUBLE,
                            attribute_space, H5P_DEFAULT);
      H5Awrite(attribute, H5T_NATIVE_DOUBLE, &max);
      H5Aclose(attribute);
      H5Sclose(attribute_space);
    }
    else{
      dataset = H5Dopen(outfile, name);
    }
    
    printf("dataset open.\n");
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d,
                                 stride_4d, count_4d, NULL);
    printf("hyperslab selected]\n");
    memspace = H5Screate_simple(rank, dimens_4d, NULL);
    printf("memspace done.\n");
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,dataspace, 
                      H5P_DEFAULT, unk);
    printf("write complete.\n");
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    
    free(unk);
  }
  return;
}
