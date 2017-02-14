/* source/physics/Eos/EosMain/Tabluated/Hdf5TableRead/eos_tabReadhdf5Tables*/

/* 
 * NAME
 *
 *  eos_tabReadhdf5
 *
 * SYNOPSIS
 *
 * void FTOC(eos_read_hdf5)(char * filename,
 *                          char * groupname,
 *                          double * temperatures,
 *                          double * densities,
 *                          double * avg_ionization,
 *                          double * ion_press,
 *                          double * ele_press,
 *                          double * ion_energy,
 *                          double * ele_energy,
 *                          double * ion_heat_capacity,
 *                          double * ele_heat_capacity,
 *                          double *entropy
 *                          //double * rosseland,
 *                          //double * planck_absorb,
 *                          //double * plank_emiss                                                                                                                                                          *                          )
 *
 *
 *
 * DESCRIPTION
 *
 *  Reads tabulated material data, eos/opacity, from an OPACPLOT datafile output (this
 *  file should be in an hdf5 file). Then writes material data to a 2 (4?) dimensional array 
 *  passed from fortran, then returns the array.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "mangle_names.h"
#include <assert.h>

void FTOC(eos_read_hdf5)(char * tableName,
			 char * groupName,
			 double * temperatures,
			 double * densities,
			 double * avg_ionization,
			 double * ion_press,
			 double * ele_press,
			 double * ion_energy,
			 double * ele_energy,
			 double * ion_heat_capacity,
			 double * ele_heat_capacity,
			 double * entropy
			 )
{
  hid_t file, dataset, dataspace, datatype, attrname, group;
  herr_t status;
  
  file = H5Fopen(tableName, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  group = H5Gopen(file, groupName);
  assert(group >= 0);

  dataset = H5Dopen(group, "temperatures");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperatures);
  dataset = H5Dopen(group, "densities");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, densities);
  dataset = H5Dopen(group, "avg_ionization");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, avg_ionization);
  dataset = H5Dopen(group, "ion_press");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ion_press);
  dataset = H5Dopen(group, "ele_press");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ele_press);
  dataset = H5Dopen(group, "ion_energy");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ion_energy);
  dataset = H5Dopen(group, "ele_energy");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ele_energy);
  dataset = H5Dopen(group, "dion_energy_dion_temp");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ion_heat_capacity);
  dataset = H5Dopen(group, "dele_energy_dele_temp");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ele_heat_capacity);

  status = H5Gclose(group);
  status = H5Fclose(file);
}
