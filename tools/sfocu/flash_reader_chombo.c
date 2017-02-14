#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5.h"
#include "flash_reader.h"
#include "options.h"

/* structs used to read bounding box dimensions from hdf5 file*/
#define PADSIZE 4
typedef struct {
	int lo_i;
	int hi_i;
	char padding[PADSIZE];
}Box1D;


typedef struct{
	int lo_i;
	int lo_j;
	int hi_i;
	int hi_j;
	char padding[PADSIZE]; /*for alignment*/
}Box2D;

typedef struct{
	int lo_i;
	int lo_j;
	int lo_k;
	int hi_i;
	int hi_j;
	int hi_k;
	char padding[PADSIZE];
}Box3D;



void err(char *err_msg){
	fprintf(stderr,"%s\n", err_msg);
	exit(1);
}



FR_Block *FR_GetBlock_HDF5_Chombo(FR_File *file, int varPos, int block_no) {
	FR_Block *out;
	char refine_level_str[MAX_STRING_LENGTH];
	hid_t handle = *( (hid_t*)(file->handle));
	int refine_level, ncells;
	int lo_i, lo_j, lo_k, hi_i, hi_j, hi_k;
	int dim, nxb, nyb, nzb;
	/*int rank;*/
	int d, b;
	/*double size;*/
	
	hid_t dspace, dset, mspace, gid ;
	herr_t status;
	hsize_t dimens, offset, count;

	dim = file->dim;
	refine_level = file->lref[block_no];
	/*size = file->size[block_no];*/
	/*rank =1;*/
	
	int lev_block0 = 0;
	for(refine_level=0; refine_level < file->lref[block_no]; refine_level++)
		lev_block0 += file->blocksPerRefineLevel[refine_level];
	b = block_no - lev_block0;
	
	sprintf(refine_level_str, "level_%d", refine_level);
	gid = H5Gopen(handle, refine_level_str);
	
	/* get offsets */
	dset = H5Dopen(gid, "data:offsets=0");
	dspace = H5Dget_space(dset);
	int *offsets = (int*)malloc((file->blocksPerRefineLevel[refine_level] + 1)*sizeof(int));
	H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, offsets);
	offset = offsets[b];
	ncells = (offsets[b+1]-offsets[b])/file->nvar;
	offset = offsets[b] + ncells*varPos;
	free(offsets);
	H5Dclose(dset);
	
	/* Number of cells in this block. Need to get coordinates of bounding box of block
	 * from the flat array into which we parsed it earlier*/
	nxb = 1; nyb = 1; nzb = 1;

	lo_i = block_no * dim * 2;
	hi_i = lo_i + dim;

	lo_j = block_no * dim * 2 + 1;
	hi_j = lo_j + dim;

	lo_k = block_no * dim * 2 + 2;
	hi_k = lo_k + dim;

	nxb = (file->bbox[hi_i]+ 1) - file->bbox[lo_i];
	if(dim > 1)
		nyb = (file->bbox[hi_j] +1) - file->bbox[lo_j];
	if(dim > 2)
		nzb = (file->bbox[hi_k] +1) - file->bbox[lo_k];
	
	/* Initialize block data struct */
	out = (FR_Block*)malloc(sizeof(FR_Block));
	double *buf = (double*)malloc(sizeof(double)*ncells);
	out->data = (double*)malloc(sizeof(double)*nxb*nyb*nzb);
	out->size[0] = nxb;
	out->size[1] = nyb;
	out->size[2] = nzb;
	int gx=0, gy=0, gz=0;
	while((nxb+2*gx)*(nyb+2*gy)*(nzb+2*gz) < ncells) {
		gx += dim > 0 ? 1 : 0;
		gy += dim > 1 ? 1 : 0;
		gz += dim > 2 ? 1 : 0;
	}
	nxb += 2*gx; nyb += 2*gy; nzb += 2*gz;
	if(nxb*nyb*nzb != ncells)
		printf("ERROR: irregular ghost cells.\n");
	
	/* Get Data */
	dset = H5Dopen(gid, "data:datatype=0");
	dspace = H5Dget_space(dset);

	dimens = ncells ; /*select all values for 1 variable for 1 block */
	count= ncells;
	mspace= H5Screate_simple(1, &dimens, NULL);
	H5Sselect_hyperslab(dspace, H5S_SELECT_SET, &offset, /*&stride*/0x0, &count, NULL);

	H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (void*)buf);
	int i, j, k;
	for(i=0; i < out->size[2]; i++)
		for(j=0; j < out->size[1]; j++)
			for(k=0; k < out->size[0]; k++)
				out->data[(i*out->size[1] + j)*out->size[0] + k] = buf[((i+gz)*nyb + j+gy)*nxb + k+gx];
	free(buf);
	H5Sclose(mspace);
	H5Sclose(dspace);
	H5Sclose(dset);
	H5Gclose(gid);
	return out;
}

FR_File* FR_open_HDF5_Chombo(hid_t handle, FR_File *out, options_t* opts){


	int nRefineLevels, numVariables, dim;
	int n_rBlocks ;
	hsize_t nB;

	hid_t gid, dspace, dset, attrib, varNameType;
	herr_t status;

	char component_name[MAX_STRING_LENGTH];
	char level_name[MAX_STRING_LENGTH];

	char tmpVarName[FR_VAR_STRING_SIZE + 1];
	

	hid_t tmpType; 

	float dx;
	
	int *buf;
	int idx_bbox; 
	
	int i, j, k;


	/*prepare file struct*/
	out->bbox = (double *)malloc(sizeof(double) );
	out->nodetype=NULL;

	/*We are not going to use these. But sfocu code checks for them so make them 
	 * trivially zero*/
	for( i=0; i < FR_MDIM; i++)
		out->ncells_vec[i] =0;


	/*Currently not doing particles with Chombo*/
	out->totalparticles=0;
	out->numRealPartProps=0;
	out->numIntPartProps=0;	
	

	out->format= FR_HDF5_CHOMBO;
	
	
	/*get total number of refinement level*/
	attrib= H5Aopen_name( handle, "num_levels");
	status= H5Aread(attrib, H5T_NATIVE_INT, (void *)&nRefineLevels);
	H5Aclose(attrib); 
	out->nLevels= nRefineLevels;
	


	/*Dimensions*/
	gid= H5Gopen(handle, "Chombo_global");
	attrib= H5Aopen_name( gid, "SpaceDim");
	status= H5Aread(attrib, H5T_NATIVE_INT, (void *)&dim);
	H5Aclose(attrib); 
	out->dim=dim;
	H5Gclose(gid);
	



	/*number of components*/
	attrib= H5Aopen_name( handle, "num_components");
	status= H5Aread(attrib, H5T_NATIVE_INT, (void *)&numVariables);
	H5Aclose(attrib); 
	out->nvar= numVariables;

	/*Get component names*/
	varNameType = H5Tcopy(H5T_C_S1);
	H5Tset_size(varNameType, 4);
	
	tmpVarName[FR_VAR_STRING_SIZE ] ='\0';


	

	for(i=0; i < numVariables; i++){
		sprintf(component_name,"component_%d", i); 
	
		attrib= H5Aopen_name( handle, component_name);
		status= H5Aread(attrib, varNameType, (void *)tmpVarName);
		H5Aclose(attrib); 


		strncpy(out->varnames[i], (char *)tmpVarName, FR_VAR_STRING_SIZE);
		out->varnames[i][FR_VAR_STRING_SIZE] = '\0';
		out->vartypes[i]= UNK;
	}

	H5Tclose(varNameType);
	

	/* Get Total Number of Blocks, refine level for each Block */
	
	out->nblocks = 0;
	out->blocksPerRefineLevel = (int *)malloc(sizeof(int) * nRefineLevels);

	

	idx_bbox=0;

	Box1D* b1=NULL;
	Box2D* b2=NULL;
	Box3D* b3=NULL;
		
	for(i =0; i < nRefineLevels; i++){
		sprintf(level_name, "level_%d", i);
		gid= H5Gopen(handle, level_name);

		dset= H5Dopen(gid, "boxes");
		dspace = H5Dget_space(dset);

		status = H5Sget_simple_extent_dims(dspace, &nB, NULL);
		n_rBlocks = (int)nB;
		
		int block0 = out->nblocks;
		out->nblocks += n_rBlocks; 
		out->blocksPerRefineLevel[i] = n_rBlocks;


		/*get bounding box for all the blocks on this refine level */
		
		/*
		attrib= H5Aopen_name( gid, "dx");
		status= H5Aread(attrib, H5T_NATIVE_FLOAT, (void *)&dx);
		H5Aclose(attrib); 
		*/

		if(dim == 1){
			b1 = (Box1D *)malloc(sizeof(Box1D) * n_rBlocks);
			buf = (int *)b1;
		}else if(dim == 2){
			b2 = (Box2D *)malloc(sizeof(Box2D) * n_rBlocks);
			buf = (int *)b2;
		}else{
			b3 = (Box3D *)malloc(sizeof(Box3D) * n_rBlocks);
			buf = (int *)b3;
		}
		tmpType= H5Dget_type(dset);
		status= H5Dread(dset, tmpType, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
		H5Tclose(tmpType);

		double *bptr = realloc(out->bbox, sizeof(double) * (out->nblocks * 2 * dim));
		if(bptr != NULL)
			out->bbox = bptr;
		else
			err("Memory Allocation failed for bbox array");
		
		/*unroll struct into the out->bbox double array*/
		for(j = 0; j < n_rBlocks; j++){
			for(k=0; k < (dim *2); k++){
				out->bbox[idx_bbox++] = (*buf);
				buf++;
			}
			buf += (PADSIZE/sizeof(int)); 
		}
			
		buf = NULL;



		/***clean up************/
		free(b1);
		free(b2);
		free(b3);
		
		H5Sclose(dspace);
		H5Dclose(dset);
		H5Gclose(gid);

	}

	/*All Chombo blocks are LEAF blocks*/
	out->nodetype = (int *) malloc(sizeof (int) * out->nblocks);
	for(j=0; j < out->nblocks; j++)
		out->nodetype[j]= FR_LEAF_NODE; 
	
	/* Fill size and coords array. This data set is redundant so we will fill it with zeros so that it trivially
	 * passes some checks in the sameblock() function from the old sfocu code for
	 * compatibility
	 */
	out->coord = (double *) malloc(sizeof(double) * out->nblocks * out->dim);
	out->size = (double *) malloc(sizeof(double) * out->nblocks * out->dim);
	for(j=0; j < (out->nblocks * out->dim ); j++){
		out->coord[j] = 0.;
		out->size[j] = 0.;
	}

	/* Fill lref array */
	int idx_lref =0;
	out->lref = (int *) malloc( sizeof(int) * out->nblocks);
	for(j = 0; j < out->nLevels; j++){
		for(k=0; k < out->blocksPerRefineLevel[j]; k++){
			out->lref[idx_lref++] = j;
		}
	}



	return out;

	

	
}

