/*****************************************************************************
* 
* This code is based on public domain code "isovis" by Mike Krogh, NCSA, 1990
* 
* He asked, but did not require, that the following message be included in all
* derived works:
* 
* Portions developed at the National Center for Supercomputing Applications at
* the University of Illinois at Urbana-Champaign.
* 
* THE UNIVERSITY OF ILLINOIS GIVES NO WARRANTY, EXPRESSED OR IMPLIED, FOR THE
* SOFTWARE AND/OR DOCUMENTATION PROVIDED, INCLUDING, WITHOUT LIMITATION,
* WARRANTY OF MERCHANTABILITY AND WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE
* 
****************************************************************************/

/*
 * This program implements the marching cubes surface tiler described by
 * Lorensen & Cline in the Siggraph 87 Conference Proceedings.
 *
 * Marching cubes portion written by Mike Krogh, NCSA, Feb.  2, 1990
 * Adapted for surface area measurement by Dean Townsley, 2011
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mangle_names.h"

/* #define DEBUG */

typedef struct {
    int    nverts;
    int    verts[8];
    int    nedges;
    int    edges[12];
    int    npolys;
    int    polys[30];
    } CELL_ENTRY;
#include "cell_table.h"

typedef struct {
    int numpolys;
    double verts[80][3][3];  /* enough for 10 polys per cell split 3 times */
    } Poly_list;


int Driver_abortFlashC(char* message);

/* prototypes */
void add_polygon(double *vert1, double *vert2, double *vert3, Poly_list *polydata );
double sum_surface_area(Poly_list* polydata,double dx,double dy,double dz);
void get_cell_verts(int index,int x1,int y1,int z1,double xtrans,double ytrans,double ztrans,
                    double threshold,double crossings[13][3]);
void get_cell_polys(int index,int *npolys,double crossings[13][3],Poly_list* polydata);
void chop_surface_cell(int x,int y,int z,int xsize,int ysize,int zsize,Poly_list* polydata);
void calc_index_and_temps(register double *data,int x1,int y1,int z1,int xdim,int ydim,int zdim,
                          register double threshold,int *index);
void print_polys(FILE* fs, Poly_list *polydata);

/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/

double DATA1,DATA2,DATA3,DATA4,DATA5,DATA6,DATA7,DATA8;
int XDIMYDIM;





/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/

/* This is the main routine.  It computes the area of the isosurface at the
   given level "threshold", where the surface is confined to the "interior"
   of the block of cells.  The interior is taken to extend out to halfway
   between the outermost two layers of cells.  This allows contiguous
   blocks of cell-centered data to have their areas evaluated without
   overlap.

   When computing surface polygons via marching cubes, the cell-centered
   data is treated like point data (vertex data) for interpolation purposes
   and then the resulting surfaces are clipped on the planes defining the
   outer faces of the interior region before adding their areas to the total.
*/

void FTOC(interior_iso_surface_area)(data,xdim,ydim,zdim,dx,dy,dz,threshold,area)
double *data;
int *xdim,*ydim,*zdim;
double *dx,*dy,*dz;
double *threshold;
double *area;
{

  register int x,y,z;
  register int xdim1,ydim1,zdim1;
  int index;
  int npolys;
  int old_verts;
  int edge[13];
  double crossings[13][3];
  Poly_list *polydata;

  polydata = malloc(sizeof(Poly_list));

  zdim1 = *zdim-1;
  ydim1 = *ydim-1;
  xdim1 = *xdim-1;

  XDIMYDIM = *xdim * *ydim;

  npolys=0;     /* keep count of total polygons */
  *area=0;

  for (z=0;z<zdim1;z++)      /* process each cell in the volume */
    for (y=0;y<ydim1;y++)
      for (x=0;x<xdim1;x++) {
#ifdef DEBUG
        fprintf(stderr, "  cell %i %i %i\n", x, y, z);
#endif
        polydata->numpolys=0;
        calc_index_and_temps(data,x,y,z,*xdim,*ydim,*zdim,*threshold,&index);
#ifdef DEBUG
        fprintf(stderr, "  temps %e %e %e %e %e %e %e %e\n",
	                DATA1,DATA2,
	                DATA3,DATA4,
	                DATA5,DATA6,
	                DATA7,DATA8);
#endif
        if (index) {
          get_cell_verts(index,x,y,z,0.0,0.0,0.0,*threshold,crossings);
          get_cell_polys(index,&npolys,crossings,polydata);
#ifdef DEBUG
        fprintf(stderr, "got %i polys before chop\n", polydata->numpolys);
	print_polys(stderr,polydata);
#endif
	  chop_surface_cell(x,y,z,*xdim,*ydim,*zdim,polydata);
#ifdef DEBUG
        fprintf(stderr, "now %i polys after chop\n", polydata->numpolys);
	print_polys(stderr,polydata);
#endif
	  *area += sum_surface_area(polydata,*dx,*dy,*dz);
        }
#ifdef DEBUG
        fprintf(stderr, "  finish, cumulative area %e\n\n", *area);
#endif
      }

  free(polydata);

}



/**************************** sum_sorface_area ********************************/

double sum_surface_area(polydata,dx,dy,dz)
Poly_list *polydata;
double dx,dy,dz;
/* sum the surface area of all polygons in the polydata list */
{
  int poly,d;
  double area;
  double deltas[3], vect1[3], vect2[3], cross[3];

  deltas[0]=dx; deltas[1]=dy; deltas[2]=dz;
  area = 0;
  for (poly=0;poly<(polydata->numpolys);poly++) {

    for (d=0;d<3;d++) {
      vect1[d] = ( (polydata->verts)[poly][1][d] - (polydata->verts)[poly][0][d] )*deltas[d];
      vect2[d] = ( (polydata->verts)[poly][2][d] - (polydata->verts)[poly][0][d] )*deltas[d];
    }

    cross[0] = vect1[1]*vect2[2]-vect1[2]*vect2[1];
    cross[1] = vect1[0]*vect2[2]-vect1[2]*vect2[0];
    cross[2] = vect1[0]*vect2[1]-vect1[1]*vect2[0];

    area += 0.5*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

  }

  return area;

}



/**************************** chop_surface_cell *******************************/

void chop_surface_cell(x,y,z,xsize,ysize,zsize,polydata)
int x,y,z,xsize,ysize,zsize;
Poly_list *polydata;
/* split polygons for face, edge or corner cells and leave only sections
   interior to halfway through the outer grid intervals */
{
  int coord[3],size[3], ndropvert, dropvert[3], nkeepvert, keepvert[3], otherdims[2];
  double cutat, intx1, intx2, inty1, inty2, newvert[2][3], newy;
  int i, dropi, odi, dim, poly, newi, polystocheck, appendedpolys;

  /* just return if we are not on a face/edge/corner */
  if ( (x>0) && (x<xsize-2) && (y>0) && (y<ysize-2) && (z>0) && (z<zsize-2) ) return;

  coord[0]=x; coord[1]=y; coord[2]=z;
  size[0]=xsize; size[1]=ysize; size[2]=zsize;
  for (dim=0;dim<3;dim++) { /* loop over dimensions */
    polystocheck=polydata->numpolys; /* we may add or remove polygons */
    poly=0;
    appendedpolys=0;
    while (poly < polystocheck ) {

      /* first determine which verticies need to be droped and where we are cutting */
      if (coord[dim]==0) {  /* lower face */
	cutat=0.5;
        ndropvert=0; nkeepvert=0;
        for (i=0;i<3;i++) {
          if ( (polydata->verts)[poly][i][dim] < 0.5 ) { dropvert[ndropvert++]=i; }
	  else { keepvert[nkeepvert++]=i; }
        }
      } else if (coord[dim]==size[dim]-2 ) { /* upper face */
        cutat = (double)size[dim]-1.5;
        ndropvert=0; nkeepvert=0;
        for (i=0;i<3;i++) {
          if ( (polydata->verts)[poly][i][dim] > cutat ) { dropvert[ndropvert++]=i; }
	  else { keepvert[nkeepvert++]=i; }
        }
      } else {
        /* break out of poly loop because not at face in this dimension */
        break;
      }

      /* now actually trim polys */
      /* other dims to interpolate in */
      if (dim==0) { otherdims[0]=1; otherdims[1]=2; }
      else if (dim==1) { otherdims[0]=0; otherdims[1]=2; }
      else if (dim==2) { otherdims[0]=0; otherdims[1]=1; }
      switch (ndropvert) {
        case 0 : /* do nothing */
	  poly++;
	  break;

	case 1 :
	  /* trim off one vertex at cut plane */
	  /* produces an extra polygon */
	  /* x here is interpolation coordinate not necessarily x dimension*/
	  /* point 1 is vertex that is being removed */
	  intx1 = (polydata->verts)[poly][dropvert[0]][dim];
	  for (newi=0;newi<2;newi++) { /* a new vertex on each cut edge */
	    intx2 = (polydata->verts)[poly][keepvert[newi]][dim];
	    /* y here is not literally y, but dependent interpolation coordinate */
	    for (odi=0;odi<2;odi++) { /* for each other dimension */
	      inty1 = (polydata->verts)[poly][dropvert[0]][otherdims[odi]];
	      inty2 = (polydata->verts)[poly][keepvert[newi]][otherdims[odi]];
	      newy = inty1+ (inty2-inty1)/(intx2-intx1)*(cutat-intx1);
	      newvert[newi][otherdims[odi]]=newy;
	    }
	    newvert[newi][dim] = cutat;
	  }
	  /* now stitch the new ond old vertices together into */
	  /* in place of current poly, put poly that has both old vertices */
	  /* just requires changing out the cut off vertex */
          memcpy( &(polydata->verts)[poly][dropvert[0]][0], newvert[0], 3*sizeof(double) );
	  /* add a poly defined by second kept vertex and new vertices */
	  memcpy( &(polydata->verts)[polydata->numpolys][0][0],
	          &(polydata->verts)[poly][keepvert[1]][0], 3*sizeof(double) );
	  memcpy( &(polydata->verts)[polydata->numpolys][1][0], newvert[0], 3*sizeof(double) );
	  memcpy( &(polydata->verts)[polydata->numpolys][2][0], newvert[1], 3*sizeof(double) );
	  (polydata->numpolys)++;

	  poly++;
	  appendedpolys++;
	  break;

        case 2 :
          /* trim two vertices by interpolating and moving them to the clip plane */
	  /* move each vertex in this coord*/
	  /* x here is interpolation coordinate not necessarily x dimension*/
	  /* point 1 is vertex that is not moving */
	  intx1 = (polydata->verts)[poly][keepvert[0]][dim];

	  for (dropi=0;dropi<2;dropi++) { /* for each vertex being dropped (moved) */
	    intx2 = (polydata->verts)[poly][dropvert[dropi]][dim];
	    /* y here is not literally y, but dependent interpolation coordinate */
	    for (odi=0;odi<2;odi++) { /* for each other dimension */
	      inty1 = (polydata->verts)[poly][keepvert[0]][otherdims[odi]];
	      inty2 = (polydata->verts)[poly][dropvert[dropi]][otherdims[odi]];
	      newy = inty1+ (inty2-inty1)/(intx2-intx1)*(cutat-intx1);
	      (polydata->verts)[poly][dropvert[dropi]][otherdims[odi]]=newy;
	    }
	    (polydata->verts)[poly][dropvert[dropi]][dim] = cutat;
	  }

	  poly++;
	  break;
	
	case 3 :
	  /* remove the polygon */
	  /* move one from the end of the list into its place */
	  if (poly==polydata->numpolys-1) {
            (polydata->numpolys)--;
	    polystocheck--;
	  } else {
	    memcpy(&(polydata->verts)[poly][0][0],
	           &(polydata->verts)[polydata->numpolys-1][0][0], 9*sizeof(double));
            (polydata->numpolys)--;
	    /* don't need to check polys that have been appended in slicing operations */
	    /* go on to next poly */
	    /* but need to account for the fact that they moved from end */
	    if (appendedpolys>0) { poly++; appendedpolys--; }
	    /* otherwise we actually dropped one, so don't check ones that aren't there */
	    /* don't increment poly index, need to check poly just moved. */
	    else { polystocheck--; }
	  }
	  break;
	
      } /* end switch on num vertices being dropped */
        

    } /* end loop over polys */
  } /* end loop over dimensions */

}


/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/

void calc_index_and_temps(data,x1,y1,z1,xdim,ydim,zdim,threshold,index)
register double *data;
int x1,y1,z1;
int xdim,ydim,zdim;
register double threshold;
int *index;
/* This subroutine calculates the index and creates some global */
/* temporary variables (for speed). */
{

  register double *tmp;

  *index = 0;

  tmp = data + (z1*XDIMYDIM) + (y1*xdim) + x1;

  *index += (threshold <= (DATA1 = *(tmp)));
  *index += (threshold <= (DATA2 = *(tmp + 1))) * 2;

  tmp += xdim;
  *index += (threshold <= (DATA3 = *(tmp + 1))) * 4;
  *index += (threshold <= (DATA4 = *(tmp))) * 8;

  tmp = tmp - xdim + XDIMYDIM;
  *index += (threshold <= (DATA5 = *(tmp))) * 16;
  *index += (threshold <= (DATA6 = *(tmp + 1))) * 32;

  tmp += xdim;
  *index += (threshold <= (DATA7 = *(tmp + 1))) * 64;
  *index += (threshold <= (DATA8 = *(tmp))) * 128;
 
}




/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/

void get_cell_verts(index,x1,y1,z1,xtrans,ytrans,ztrans,threshold,crossings)
int index;
int x1,y1,z1;
double xtrans,ytrans,ztrans;
double threshold;
double crossings[13][3];
{

  register int i;
  register int x2,y2,z2;
  int nedges;
  int crnt_edge;
 
#define linterp(a1,a2,a,b1,b2) ((double)(((a-a1) * (double)(b2-b1) / (a2-a1)) + (double)b1))

  x2 = x1+1;
  y2 = y1+1;
  z2 = z1+1;

  nedges = cell_table[index].nedges;
  for (i=0;i<nedges;i++) {
     crnt_edge = cell_table[index].edges[i];
     switch (crnt_edge) {
	case 1:
		crossings[1][0] = linterp(DATA1,DATA2,threshold,x1,x2)+xtrans;
		crossings[1][1] = (double)y1+ytrans;
		crossings[1][2] = (double)z1+ztrans;
		break;

	case 2:
		crossings[2][1] = linterp(DATA2,DATA3,threshold,y1,y2)+ytrans;
		crossings[2][0] = (double)x2+xtrans;
		crossings[2][2] = (double)z1+ztrans;
		break;

	case 3:
		crossings[3][0] = linterp(DATA4,DATA3,threshold,x1,x2)+xtrans;
		crossings[3][1] = (double)y2+ytrans;
		crossings[3][2] = (double)z1+ztrans;
		break;

	case 4:
		crossings[4][1] = linterp(DATA1,DATA4,threshold,y1,y2)+ytrans;
		crossings[4][0] = (double)x1+xtrans;
		crossings[4][2] = (double)z1+ztrans;
		break;

	case 5:
		crossings[5][0] = linterp(DATA5,DATA6,threshold,x1,x2)+xtrans;
		crossings[5][1] = (double)y1+ytrans;
		crossings[5][2] = (double)z2+ztrans;
		break;

	case 6:
		crossings[6][1] = linterp(DATA6,DATA7,threshold,y1,y2)+ytrans;
		crossings[6][0] = (double)x2+xtrans;
		crossings[6][2] = (double)z2+ztrans;
		break;

	case 7:
		crossings[7][0] = linterp(DATA8,DATA7,threshold,x1,x2)+xtrans;
		crossings[7][1] = (double)y2+ytrans;
		crossings[7][2] = (double)z2+ztrans;
		break;

	case 8:
		crossings[8][1] = linterp(DATA5,DATA8,threshold,y1,y2)+ytrans;
		crossings[8][0] = (double)x1+xtrans;
		crossings[8][2] = (double)z2+ztrans;
		break;

	case 9:
		crossings[9][2] = linterp(DATA1,DATA5,threshold,z1,z2)+ztrans;
		crossings[9][1] = (double)y1+ytrans;
		crossings[9][0] = (double)x1+xtrans;
		break;

	case 10:
		crossings[10][2] = linterp(DATA2,DATA6,threshold,z1,z2)+ztrans;
		crossings[10][1] = (double)y1+ytrans;
		crossings[10][0] = (double)x2+xtrans;
		break;

	case 11:
		crossings[11][2] = linterp(DATA4,DATA8,threshold,z1,z2)+ztrans;
		crossings[11][1] = (double)y2+ytrans;
		crossings[11][0] = (double)x1+xtrans;
		break;

	case 12:
		crossings[12][2] = linterp(DATA3,DATA7,threshold,z1,z2)+ztrans;
		crossings[12][1] = (double)y2+ytrans;
		crossings[12][0] = (double)x2+xtrans;
		break;

    } /* end switch */
  } /* end for */
}




/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/

void get_cell_polys(index,npolys,crossings,polydata)
int index;
int *npolys;
double crossings[13][3];
Poly_list * polydata;
/* This subroutine will calculate the polygons */
{

  register int num_o_polys;
  register int poly;
  double *p1,*p2,*p3;
  double n1[3],n2[3],n3[3];

  num_o_polys = cell_table[index].npolys;
  for (poly=0;poly<num_o_polys;poly++) {

    p1 = &crossings[cell_table[index].polys[(poly*3)]][0];
    p2 = &crossings[cell_table[index].polys[(poly*3)+1]][0];
    p3 = &crossings[cell_table[index].polys[(poly*3)+2]][0];

    add_polygon(p1,p2,p3,polydata);

  }

  (*npolys) += num_o_polys;

}


void add_polygon(double *vert1, double *vert2, double *vert3, Poly_list *polydata )
{

  if (polydata->numpolys==80) {
    Driver_abortFlashC("Not enough space for polygons for surface area\n");
  }
  memcpy( &(polydata->verts)[polydata->numpolys][0][0], vert1, 3*sizeof(double));
  memcpy( &(polydata->verts)[polydata->numpolys][1][0], vert2, 3*sizeof(double));
  memcpy( &(polydata->verts)[polydata->numpolys][2][0], vert3, 3*sizeof(double));

  polydata->numpolys++;
}


void print_polys(FILE* fs, Poly_list *polydata)
{
  int poly, v, d;
  if (polydata->numpolys==0) {
    fprintf(fs,"no polygons\n");
    return;
  }

  for (poly=0;poly<polydata->numpolys;poly++) {
    fprintf(fs,"polygon %i:\n", poly);
    for (v=0;v<3;v++) {
      fprintf(fs," vertex %i:",v);
      for (d=0;d<3;d++) {
        fprintf(stderr," %e", (polydata->verts)[poly][v][d] );
      }
      fprintf(fs,"\n");
    }
  }

  return;
}
