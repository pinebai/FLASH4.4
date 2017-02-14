/*

NAME

  eos_tabBrowsehdf5Tables.c

SYNOPSIS

  call eos_tabBrowsehdf5C(char * tableName,
                                    char * groupname,
				    int * ntemp,
				    int * ndens,
				    int * ngroups)

DESCRIPTION

  This program contains a function to be called from eos_tabReadOpacplot.F90. It first checks whether the 
  specified file is indeed an hdf5 file, and if true, proceeds to open the file. The function, 
  eostabBrowseOpacplotTablesC, then opens a specified group, collects attribute data (ntemp, ndens,
  ngroups), then closes the group and file. 

ARGUMENTS

NOTES


*/


#include <stdio.h>
#include <stdlib.h>
#include "mangle_names.h"
#include "hdf5.h"
#include <assert.h>

void FTOC(eos_tabbrowsehdf5)(char * tableName,
			     char * groupName,
			     int * ntemp,
			     int * ndens)
{

  hid_t attrname, oppf, group;
  herr_t status;

  printf("Table name = %s, Group name = %s\n", tableName, groupName);

  oppf = H5Fopen(tableName, H5F_ACC_RDONLY, H5P_DEFAULT);   //open file, tableName
  assert (oppf >= 0);

  group = H5Gopen(oppf, groupName);          //open groupName
  assert (group >= 0);
  
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8)
    /* This function is now deprecated by HDF5 */
#  define H5Aopen(ds,name,ignore) H5Aopen_name(ds,name)
#endif

  attrname = H5Aopen(group, "ntemp", H5P_DEFAULT);         // open attribute ntemp
  status = H5Aread(attrname, H5T_NATIVE_INT, ntemp);        // read attribute ntemp
  
  attrname = H5Aopen(group, "ndens", H5P_DEFAULT);   // open attribute ndens
  status = H5Aread(attrname, H5T_NATIVE_INT, ndens);    // read attribute ndens

  //  attrname = H5Aopen(group, "ngroups", H5P_DEFAULT);   // open attribute ngroups
  //  status = H5Aread(attrname, H5T_NATIVE_INT, ngroups); // read attribute ngroups
  
  status = H5Fclose(oppf);
  status = H5Gclose(group);

}
