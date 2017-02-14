#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <ctype.h>
#include "flash_reader.h"
#include "sameblock.h"
#include "namecmp.h"
#include "options.h"

#define SQUARE(X) (X)*(X)
/* FIXME:
   
   When the meshes have no blocks in common, we read unintialized data when
   we print the report

   TODO: would be nice to have histograms (bin errors over cells)
*/

#define NO_BLOCK -1 /* Not one of paramesh's nodetypes, I think */


/* norm(a, b)

   Same norm as FOCU; if you change this change the output message below as
   well. This is called only when a!=b 
*/

double norm(double a, double b){
  return fabs(2*(a-b))/(fabs(a+b) > 1e-99 ? fabs(a+b) : 1e-99); 
}

/* l2norm(a,b)

   Computes an L2-norm.  Should only work when a Bad Block detected.

   side effects a difference sum, we can fetch that later.
 */
/*
double lx_errorNorm(double a, double b, double order){
  
  double error = 0.;
  
  if (order != 0.0)
    error = pow( pow(abs(b-a), order), 1./order);
  
  return error;
}
*/
/*
double lx_relativeNorm(double a, double b, double order){
  
  double error = 0.;

  if(order != 0.0 && a != 0.0)
    error = pow( (pow(abs(b - a), order)) / pow(abs(a), order)), 1./order);

  return error;
}
*/


char* strip(char* s) {
  /**
   * strip all whitespace (leading, trailing, and internal)
   * from string 's' and return the stripped string
   */
  int idx, idx2, charCount;
  char* newStr;

  /* tally up all non-whitespace characters in 's' */
  idx = 0;
  charCount = 1; /* 1 automatic character for the null-terminator */
  while (s[idx] != 0) {
    if (!isspace(s[idx])) {
      charCount++;
    }
    idx++;
  }

  /* malloc memory for all non-whitespace characters in 's' */
  newStr = (char *) malloc(sizeof(char) * charCount);
  idx  = 0;
  idx2 = 0;
  while (s[idx] != 0) {
    if (!isspace(s[idx])) {
      newStr[idx2] = tolower(s[idx]);
      idx2++;
    }
    idx++;
  }
  /* add the null terminator */
  newStr[idx2] = 0;
  return newStr;
}

static double *_pt_qsort_partbuf;
static int _pt_qsort_numprops, _pt_qsort_proptag;

static int _pt_qsort_compare(const void *_ai, const void *_bi) {
  int ai = *(const int*)_ai;
  int bi = *(const int*)_bi;
  double a = _pt_qsort_partbuf[_pt_qsort_proptag + _pt_qsort_numprops*ai];
  double b = _pt_qsort_partbuf[_pt_qsort_proptag + _pt_qsort_numprops*bi];
  return a < b ? -1 : b < a ? 1 : 0;
}

/* sorts particle indices by tag value, returns indices array */
int* sort_particle_tags(double *partbuf, int numparts, int numprops, int proptag) {
  int *perm = (int*)malloc(numparts*sizeof(int));
  int i;
  for(i=0; i < numparts; i++)
    perm[i] = i;
  _pt_qsort_partbuf = partbuf;
  _pt_qsort_numprops = numprops;
  _pt_qsort_proptag = proptag;
  qsort(perm, numparts, sizeof(int), _pt_qsort_compare);
  return perm;
}

/* returns particle index in partbuf of particle with given tag */
int find_particle_tag(double tag, double *partbuf, int *sortperm, int numparts, int numprops, int proptag) {
  int a = 0, b = numparts;
  while(a < b) {
    int m = a + (b-a)/2;
    double tag1 = partbuf[proptag + sortperm[m]*numprops];
    if(tag < tag1)
      b = m;
    else if(tag > tag1)
      a = m + 1;
    else
      return sortperm[m];
  }
  return -1;
}

/* sfocu(argc, argv)

   Compares two flash checkpoint files, to decide wether or not they have
   the same data. Does not load all the file into memory at once, so should
   be able to deal with large datasets.
*/
int sfocu(int argc, char **argv){
  FR_File *A, *B, *C;
  FR_Block *blA, *blB;
  FR_ParticlesAllProps *pAllA, *pAllB;
  int i, blockA, blockB, partIndexA, partIndexB, tagIndexA, tagIndexB, dim, failed, firstTime, a, b;
  int blkIndexA, blkIndexB;
  int particlesRead, particlesGotten;
  int *AtoB; /* AtoB[blockA]==blockB */
  int *pAtoB; /* pAtoB[particleA]==particleB */

  double error = 0.0;
  int ivar;
  int iprop;
  char normChar ='1';
  double propValueA, propValueB;
  int cellCount; /*a count to use for iterating over different grid variable types*/

  struct{
    int blocksCompared; /* number of blocks we compare */
    int extraBlocksA; /* leaf blocks in A but not B */
    int extraBlocksB;
    int nvar; /* number of variables we compare */
    int fail; /* if false, the files are "identical" */
    int badBlocks[FR_MAXVARS]; /* num leaf blocks where a var differs by > err_tol */
    double sumA[FR_MAXVARS], sumB[FR_MAXVARS];
    double maxA[FR_MAXVARS], maxB[FR_MAXVARS];
    double minA[FR_MAXVARS], minB[FR_MAXVARS];
    double maxError[FR_MAXVARS]; /* max of error norm */
    double maxErrorVal[2][FR_MAXVARS]; 
    double minError[FR_MAXVARS]; 
    double minErrorVal[2][FR_MAXVARS];  /*we calculate but don't report this */
    double absError[FR_MAXVARS]; /* sup |a-b| */
    double absErrorVal[2][FR_MAXVARS];
    double magError[FR_MAXVARS]; /* sup |a-b| / max( sup |a|, sup |b|, 1e-99) */
    double L2Norm[FR_MAXVARS]; /*maximum L2-Norm calculation for comparison */
    double sumABDiff[FR_MAXVARS]; /*sum of differences used for L2-calculation*/
    double sumSqVar[FR_MAXVARS]; /*sum of values for Ln-Error norm calculation*/
    char varnames[FR_MAXVARS][FR_VAR_STRING_SIZE+1]; /* variables we compare */
  } report = {0,0,0,0,0};

  struct{
    int particlesCompared; /* number of particles we compare */
    int extraParticlesA; /* particles in A but not B */
    int extraParticlesB;
    int numRealPartProps; /* number of properties we compare */
    int fail; /* if false, the files are "identical" */
    int badParticles[FR_MAXVARS]; /* number of particles where a prop differs by > perr_tol */
    double sumA[FR_MAXVARS], sumB[FR_MAXVARS];
    double maxA[FR_MAXVARS], maxB[FR_MAXVARS];
    double minA[FR_MAXVARS], minB[FR_MAXVARS];
    double maxError[FR_MAXVARS]; /* max of error norm */
    double maxErrorVal[2][FR_MAXVARS]; 
    double minError[FR_MAXVARS]; 
    double minErrorVal[2][FR_MAXVARS];  /*we calculate but don't report this*/
    double absError[FR_MAXVARS]; /* sup |a-b| */
    double absErrorVal[2][FR_MAXVARS];
    double magError[FR_MAXVARS]; /* sup |a-b| / max( sup |a|, sup |b|, 1e-99)*/
    char realPartPropNames[FR_MAXVARS][FR_PART_PROP_STRING_SIZE+1]; /* properties we compare */
    int ipropA[FR_MAXVARS], ipropB[FR_MAXVARS]; /* maps from report props to actual file props */
  } particleReport = {0,0,0,0,0};
  
  options_t opts;

  if(parse_cmdline(&opts, argc, argv)){
    /*printf("Usage: %s [-hvi] [-t mesh_tol] [-e err_tol] file1 file2\n", argv[0]); */
    printf("Usage: %s [-hv] [-t mesh_tol] [-e err_tol] [-p perr_tol] [-I var1 [-I var2 [..]]] [-s svar1 [-s svar2 [..]]] file1 file2\n", argv[0]);
    printf("  -h: This help screen\n");
    printf("  -v: Verbose output\n");
    printf("  -Ln: Show Ln error norms, default is -L1\n");
    printf("  -bA or -bB: Consider file A or B the benchmark\n");
    printf("  -r: Reorder given files if numbers of blocks are different\n");
    printf("  -e: Differences in mag-error value of less than err_tol\n");
    printf("      will not cause a failure. Default value is 0.0\n");
    printf("  -p: Differences in particle mag-error value of less than perr_tol\n");
    printf("      will not cause a failure. Default value is 0.0\n");
    /*    printf("  -i: Ignore differences in refinement levels between corresponding blocks,\n");
          printf("      continuing with comparison as long as blocks contain same # of zones\n");*/
    printf("  -t: Set mesh tolerance\n");
    printf("      Default value is 0.0\n");
    printf("  -I: Specify name(s) of variable(s) to ignore, max length 4 characters each - use '*' for all.\n");
    printf("      These options are processed BEFORE any -s options.\n");
    printf("  -s: Specify additional name(s) of variable(s) to compare, max length 4 characters each.\n");
    printf("      Use this to force recognition of SCRATCH vars not kept in the FLASH3 'UNK' array.\n");
    printf("  -a: Sfocu will atttempt to locate all grid-based variables from HDF5's metadata instead\n");
    printf("      of 'unknown names'.  This should allow all variables in a file to be read, regardless\n");
    printf("      of type.  This feature is still experimental.\n");
    return 1;
  }
  
  /* Open files */
  if(!(A = FR_open(opts.file1, opts))){
      printf("Couldn't read file: %s\n", opts.file1);
      return 1;
  } 

  if(!(B = FR_open(opts.file2, opts))){
    printf("Couldn't read file: %s\n", opts.file2);
    return 1;
  }
  
  if(opts.err_tol < 0) {
    printf("Mag-error tolerance may not be less than zero\n");
    return 1;
  }

  if(opts.perr_tol < 0) {
    printf("Particle mag-error tolerance may not be less than zero\n");
    return 1;
  }

  switch(opts.norm_order) {
  case (1):
    normChar = '1';
    break;
  case (2):
    normChar = '2';
    break;
  default:
    normChar = 'n';
    break;
  }

  if (opts.reorder) {
    int doswap = 0;
    if (opts.benchmark=='A') {
      /* want nblocksA >= nblocksB */
      if (A->nblocks < B->nblocks){
        doswap = 1;
      }
    } else {
      /* want nblocksA <= nblocksB */
      if (A->nblocks > B->nblocks){
        doswap = 1;
      }
    }
    if (doswap) {
      C=A; A=B; B=C;
    }
  }

  printf("\nA: %s\nB: %s\n\n", A->filename, B->filename);

  switch(A->format){
          case (FR_HDF5_PMESH):
                printf("File Format A: Paramesh HDF5\n");
                break;
          case (FR_HDF5_CHOMBO):
               printf("File Format A: Chombo HDF5\n");
               break;
          case (FR_HDF4):
                printf("File Format A: HDF4\n");
                break;
          case (FR_NCDF):
                printf("File Format A: NCDF\n");
                break;
  }
  switch(B->format){
          case (FR_HDF5_PMESH):
                printf("File Format B: Paramesh HDF5\n");
                break;
          case (FR_HDF5_CHOMBO):
               printf("File Format B: Chombo HDF5\n");
               break;
          case (FR_HDF4):
                printf("File Format B: HDF4\n");
                break;
          case (FR_NCDF):
                printf("File Format B: NCDF\n");
                break;
  }

  printf("Min Error: inf(2|a-b| / max(|a+b|, 1e-99) )\n");
  printf("Max Error: sup(2|a-b| / max(|a+b|, 1e-99) )\n");
  printf("Abs Error: sup|a-b|\n");
  printf("Mag Error: sup|a-b| / max(sup|a|, sup|b|, 1e-99)\n");
  printf("L%c-ErrNorm Rel/%c: ||B-A|| / max(||%c||, 1e-99)  where ||.|| is discrete L%d norm\n\n",
         normChar,opts.benchmark,opts.benchmark,opts.norm_order);


                        
          /* Check that cells per block match */
          for(i=0; i<FR_MDIM; i++) {
            if (A->ncells_vec[i]!=B->ncells_vec[i]){
              printf("Block shapes do not match: [%d,%d,%d] != [%d,%d,%d]\n",
                     A->ncells_vec[0], A->ncells_vec[1], A->ncells_vec[2],
                     B->ncells_vec[0], B->ncells_vec[1], B->ncells_vec[2]);
              printf("FAILURE\n");
              return 1;
            }
          }

        /* For Chombo Hdf5 files check that both files have similar block structure 
         * for the refine levels we will compare */
        if( A->format == FR_HDF5_CHOMBO && B->format == FR_HDF5_CHOMBO){
                

                          for(i=0; i < A->nLevels; i++)
                                if(A->blocksPerRefineLevel[i] != B->blocksPerRefineLevel[i]){
                                        printf("FAILURE\n");    
                                        return 1;
                                }   
        }
         else if (A->format == FR_HDF5_CHOMBO || B->format == FR_HDF5_CHOMBO){ 
                /*only one file is in chombo fmt*/              
                printf("Chombo checkpoint files can only be compared with another Chombo checkpoint files\n");

                return 1;
        }


        

  /* Report block shape */
  if( A->format != FR_HDF5_CHOMBO){
        printf("Block shapes for both files are: [%d,%d,%d]\n",
                A->ncells_vec[0], A->ncells_vec[1], A->ncells_vec[2]);
  }
  printf("Mag-error tolerance: %g\n", opts.err_tol);
  printf("Particle mag-error tolerance: %g\n", opts.perr_tol);

  /* Macros use this variable name */
  dim = A->dim; 

  /* Now use coordinate and size arrays to build a table translating
     leaf block numbers from one file to the other. We also use this chance
     to dismiss A's non-leaf blocks.

     No need for fancy algorithms since the mapping is close to the identity?
  */

  AtoB =(int *) malloc(sizeof(int)*A->nblocks);
  for(blockA=0; blockA<A->nblocks; blockA++){
    if (A->nodetype[blockA] != FR_LEAF_NODE){
      AtoB[blockA] = NO_BLOCK;
      continue;
    }
    
    failed = 1;
    for(i=0; i < B->nblocks; i++) {
      /* blockB test sequence = blockA + (0,-1,1,-2,2,-3,3,...) */
      blockB = (blockA + B->nblocks + ((i+1)>>1)*(i&1 ? -1 : 1)) % B->nblocks;
      if(sameblock(dim, &opts, A, B, blockA, blockB)){
        AtoB[blockA]=blockB;
        report.blocksCompared++;
        failed = 0;
        break;
      }
    }
    
    if(failed){
      AtoB[blockA] = NO_BLOCK;
      report.extraBlocksA++;
      if(opts.verbose)
        printf("No match for block %d in %s\n", blockA, A->filename);
    }
  } /* end of table-building loop */


  /* Now check for leaf blocks in B that didn't get into the translation
     table. We'll write over B's nodetype array to mark blocks we already know.
  */

  for(blockA=0; blockA<A->nblocks; blockA++)
    if (AtoB[blockA]!=NO_BLOCK)
      B->nodetype[AtoB[blockA]]=NO_BLOCK;
  
  for(blockB=0; blockB<B->nblocks; blockB++)
    if (B->nodetype[blockB]==FR_LEAF_NODE){
      report.extraBlocksB++;
      if(opts.verbose)
        printf("No match for block %d in %s\n", blockB, B->filename);
    }

  if(report.extraBlocksA || report.extraBlocksB)
    report.fail = 1;

  /* Remove variable names given in -I options from the lists generated
     from reading "unknown names" in both file A and file B - if
     they are present.  Superfluous -I options are ignored. - KW
     */
  for(i=0; i < opts.n_ignorableVarnames; i++){

    if(strcmp(strip(opts.ignorableVarnames[i]), "*") == 0) {
      A->nvar = 0;
      B->nvar = 0;
      break;
    }
    for(a=0; a < A->nvar; a++){
      if(namecmp(strip(A->varnames[a]), strip(opts.ignorableVarnames[i])) == 0) {
        for(; a < A->nvar-1; a++) {
          strncpy(A->varnames[a], A->varnames[a+1], FR_VAR_STRING_SIZE+1);
        }
        A->nvar--;
        break;
      }
    }
    for(b=0; b < B->nvar; b++){
      if(namecmp(strip(B->varnames[b]), strip(opts.ignorableVarnames[i])) == 0) {
        for(; b < B->nvar-1; b++) {
          strncpy(B->varnames[b], B->varnames[b+1], FR_VAR_STRING_SIZE+1);
        }
        B->nvar--;
        break;
      }
    }
  }

  /* Append extra variable names from -s options to the lists generated
     from reading "unknown names" in both file A and file B - but only if
     they are not present already.
     After this, the code will act exactly as if "unknown name" had listed
     the extra names.
     Note that there is no checking whether variables of the given names
     are ACTUALLY in the files - user beware. - KW
     */
  for(i=0; i < opts.n_extra_varnames; i++){
    int found = 0;
    for(a=0; a < A->nvar; a++){
      if(namecmp(strip(A->varnames[a]), strip(opts.extra_varnames[i])) == 0) {
        found = 1;
        break;
      }
    }
    if (!found)
      strcpy(A->varnames[A->nvar++], opts.extra_varnames[i]);
    found = 0;
    for(b=0; b < B->nvar; b++){
      if(namecmp(strip(B->varnames[b]), strip(opts.extra_varnames[i])) == 0) {
        found = 1;
        break;
      }
    }
    if (!found)
      strcpy(B->varnames[B->nvar++], opts.extra_varnames[i]);
  }

  /* Now check what variables are in both arrays */
  report.nvar = 0;
  for(a=0; a < A->nvar; a++){
    failed = 1;
    for(b=0; b < B->nvar; b++){
      if(namecmp(strip(A->varnames[a]), strip(B->varnames[b])) == 0){
        failed = 0;
        strcpy(report.varnames[report.nvar++], A->varnames[a]);
        break;
      }
    }
    if(failed){
      report.fail = 1;
      printf("%s present in %s but not in %s\n", A->varnames[a], A->filename, 
             B->filename);
    }
  }
  
  for(b=0; b<B->nvar; b++){
    failed = 1;
    for(a=0; a<A->nvar; a++){
      if(namecmp(strip(A->varnames[a]), strip(B->varnames[b])) == 0) {
        failed = 0;
        break;
      }
    }
    if(failed){
      report.fail = 1;
      printf("%s present in %s but not in %s\n", B->varnames[b], B->filename, 
             A->filename);
    }
  }

  /* Now check for common particle properties */
  particleReport.numRealPartProps = 0;

  /* Check that numbers of particles the same */
  if (A->totalparticles != B->totalparticles) {
    particleReport.fail = 1;
    printf("\nFiles have different number of particles:\n");
    printf("%d in %s, %d in %s\n", A->totalparticles, A->filename,
           B->totalparticles, B->filename);
  }
  else {
    if (A->totalparticles > 0 && B->totalparticles > 0) {
      for(a=0; a<A->numRealPartProps; a++) {
        failed = 1;
        for(b=0; b<B->numRealPartProps; b++){
          if(strcmp(strip(A->realPartPropNames[a]), strip(B->realPartPropNames[b])) == 0){
            failed = 0;
            particleReport.ipropA[particleReport.numRealPartProps] = a;
            particleReport.ipropB[particleReport.numRealPartProps] = b;
            strcpy(particleReport.realPartPropNames[particleReport.numRealPartProps++], A->realPartPropNames[a]);
            break;
          }
        }
        if(failed) {
          particleReport.fail = 1;
          printf("particle property '%s' present in %s but not in %s\n", A->realPartPropNames[a], A->filename, 
                 B->filename);
        }
      }

      for(b=0; b<B->numRealPartProps; b++) {
        failed = 1;
        for(a=0; a<A->numRealPartProps; a++){
          if(strcmp(strip(A->realPartPropNames[a]), strip(B->realPartPropNames[b])) == 0){
            failed = 0;
            break;
          }
        }
        if(failed){
          particleReport.fail = 1;
          printf("particle property '%s' present in %s but not in %s\n", B->realPartPropNames[b], B->filename, 
                 A->filename);
        }
      }

      pAllA = FR_GetParticlesAllProps(A, 0, A->totalparticles, &particlesGotten);
      pAllB = FR_GetParticlesAllProps(B, 0, B->totalparticles, &particlesGotten);
      pAtoB =(int *) malloc(sizeof(int)*A->totalparticles);
      
      /* find the particle tag index */
      for (tagIndexA = 0; tagIndexA < FR_MAXVARS; tagIndexA++) {
        if (strncmp(pAllA->realPropsNames[tagIndexA], "tag", 3) == 0) {
          printf("A tag index %d\n", tagIndexA);
          break;
        }
      }

      for (tagIndexB = 0; tagIndexB < FR_MAXVARS; tagIndexB++) {
        if (strncmp(pAllB->realPropsNames[tagIndexB], "tag", 3) == 0) {
          printf("B tag index %d\n", tagIndexB);
          break;
        }
      }
  
      /*find the particle blk index 
        we don't want to compare this property, because simulations on different
        number of procs will have particles with different blk properties even
        though the particles are in the same physical domain. */

      for (blkIndexA = 0; blkIndexA < FR_MAXVARS; blkIndexA++) {
        if (strncmp(pAllA->realPropsNames[blkIndexA], "blk", 3) == 0) {
          printf("A blk index %d\n", blkIndexA);
          break;
        }
      }

      for (blkIndexB = 0; blkIndexB < FR_MAXVARS; blkIndexB++) {
        if (strncmp(pAllB->realPropsNames[blkIndexB], "blk", 3) == 0) {
          printf("B blk index %d\n", blkIndexB);
          break;
        }
      }

      /* compute sorted permutation of B particle tags for fast lookup */
      int *pBsort = sort_particle_tags(
        pAllB->realProps, B->totalparticles, B->numRealPartProps, tagIndexB
      );
      
      for(partIndexA=0; partIndexA<A->totalparticles; partIndexA++) {
#if 1 /* use particle lookup by sorted tag */
        partIndexB = find_particle_tag(
          pAllA->realProps[partIndexA*pAllA->numRealProps + tagIndexA],
          pAllB->realProps, pBsort, B->totalparticles, B->numRealPartProps, tagIndexB
        );
        if(partIndexB >= 0 && 
          pAllA->realProps[partIndexA*pAllA->numRealProps + tagIndexA]
          != pAllB->realProps[partIndexB*pAllB->numRealProps + tagIndexB]) {
          printf("BUG in particle sort!\n");
        }
        pAtoB[partIndexA] = partIndexB;
        failed = partIndexB < 0;
#else /* do not use particle tag sorting */
        /* this is common */
        if (pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA] == 
            pAllB->realProps[(partIndexA * pAllB->numRealProps) + tagIndexB]) {
          pAtoB[partIndexA] = partIndexA;
          continue;
        }
        /*
        printf("partIndexA is: %d\n", partIndexA);
        printf("pAllA->numRealProps is: %d\n", pAllA->numRealProps);
        printf("tagIndexA is: %d\n", tagIndexA);
        printf("pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA] is: %e\n", pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA]);
        printf("pAllB->numRealProps is: %d\n", pAllB->numRealProps);
        printf("tagIndexB is: %d\n", tagIndexB);
        printf("pAllB->realProps[(partIndexA * pAllB->numRealProps) + tagIndexB] is: %e\n", pAllB->realProps[(partIndexA * pAllB->numRealProps) + tagIndexB]);
        */ 

        /* this is expensive but rare */
        failed = 1;
        for(partIndexB=0; partIndexB < B->totalparticles; partIndexB++) {
          if (pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA] == 
              pAllB->realProps[(partIndexB * pAllB->numRealProps) + tagIndexB]){
            pAtoB[partIndexA] = partIndexB;
            failed = 0;
            break;
          }
        }
#endif
        if(failed){
          pAtoB[partIndexA] = NO_BLOCK;
          particleReport.extraParticlesA++;
          /*      if(opts.verbose) */
          printf("No match for particle %d (tag=%16.16g) in %s\n",
            partIndexA,
            pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA],
            A->filename);
          particleReport.fail = 1;
          /*exit(1);*/
        }
      } /* end of particle table-building loop */
      
      free(pBsort);

      /* Check that particle real properties are the same */ 
      particleReport.particlesCompared = 0;
      for (i = 0; i < particlesGotten; i++) {
        if (pAtoB[i] != NO_BLOCK) {
          for(iprop=0; iprop<particleReport.numRealPartProps; iprop++) {
            if(iprop != blkIndexA) {

              propValueA = pAllA->realProps[i*A->numRealPartProps + particleReport.ipropA[iprop]];
              propValueB = pAllB->realProps[pAtoB[i]*B->numRealPartProps + particleReport.ipropB[iprop]];
              
              if (particleReport.particlesCompared == 0) {
                particleReport.sumA[iprop] = 0;
                particleReport.maxA[iprop] = propValueA;
                particleReport.minA[iprop] = propValueA;
                
                particleReport.sumB[iprop] = 0;
                particleReport.maxB[iprop] = propValueB;
                particleReport.minB[iprop] = propValueB;
                
                if(propValueA != propValueB)
                  particleReport.minError[iprop] = norm(propValueA, propValueB);
                else
                  particleReport.minError[iprop] = 0;

                particleReport.badParticles[iprop] = 0;
                particleReport.maxError[iprop] = 0.0;
                particleReport.absError[iprop] = 0.0;
              }

              if (propValueA != propValueB) {

                particleReport.badParticles[iprop]++;

                if (opts.verbose)
                  printf(" %d %s\n", i, particleReport.realPartPropNames[iprop]);
                error = norm(propValueA, propValueB);

                if(error > particleReport.maxError[iprop]) {
                  particleReport.maxError[iprop] = error;
                  particleReport.maxErrorVal[0][iprop] = propValueA;
                  particleReport.maxErrorVal[1][iprop] = propValueB;
                }
              }

              if(error < particleReport.minError[iprop]) {
                particleReport.minError[iprop] = error;
                particleReport.minErrorVal[0][iprop] = propValueA;
                particleReport.minErrorVal[1][iprop] = propValueB;
              }
              
              error = fabs(propValueA-propValueB);
              if(error > particleReport.absError[iprop]) {
                particleReport.absError[iprop] = error;
                particleReport.absErrorVal[0][iprop] = propValueA;
                particleReport.absErrorVal[1][iprop] = propValueB;
              }
              
              particleReport.sumA[iprop] += propValueA; 
              particleReport.maxA[iprop] = (propValueA > particleReport.maxA[iprop]) ? propValueA :
                particleReport.maxA[iprop];
              particleReport.minA[iprop] = (propValueA < particleReport.minA[iprop]) ? propValueA :
                particleReport.minA[iprop];
              
              particleReport.sumB[iprop] += propValueB; 
              particleReport.maxB[iprop] = (propValueB > particleReport.maxB[iprop]) ? propValueB :
                particleReport.maxB[iprop];
              particleReport.minB[iprop] = (propValueB < particleReport.minB[iprop]) ? propValueB :
                particleReport.minB[iprop];
              
            } /*if iprop != blk prop */
          }

        particleReport.particlesCompared++;
        } /* if pAtoB[i] is valid */

      }

      FR_DeleteParticlesAllProps(pAllA);
      FR_DeleteParticlesAllProps(pAllB);

    }    

    for(i=0; i<particleReport.numRealPartProps; i++) {
      error = fabs(particleReport.maxA[i]);
      error = fabs(particleReport.minA[i]) > error ? fabs(particleReport.minA[i]) : error;
      error = fabs(particleReport.maxB[i]) > error ? fabs(particleReport.maxB[i]) : error;
      error = fabs(particleReport.minB[i]) > error ? fabs(particleReport.minB[i]) : error;
      error = error > 1.e-99 ? error : 1.e-99;
      particleReport.magError[i] = particleReport.absError[i] / error;
      if (particleReport.magError[i] > opts.perr_tol) {
        particleReport.fail = 1;
      }
    }
  }

  /* Finally, check where the data doesn't match */
  for(ivar=0; ivar<report.nvar; ivar++){
    double r;
    if(opts.verbose)
      printf("\n%s's bad blocks:", report.varnames[ivar]);

    report.badBlocks[ivar] = 0;
    report.maxError[ivar] = 0.0;
    report.absError[ivar] = 0.0;
    
    firstTime = 1;


    for(blockA=0; blockA<A->nblocks; blockA++){
      if (AtoB[blockA]==NO_BLOCK)
        continue;

      if(A->format != FR_HDF5_CHOMBO){
        blA = FR_GetBlock(A, report.varnames[ivar], blockA);
        blB = FR_GetBlock(B, report.varnames[ivar], AtoB[blockA]);
      }else{
        blA = FR_GetBlock_HDF5_Chombo(A, ivar, blockA);
        blB = FR_GetBlock_HDF5_Chombo(B, ivar, AtoB[blockA]);
      }

      if(firstTime){
        firstTime = 0;

        report.sumA[ivar] = 0;
        report.maxA[ivar] = blA->data[0];
        report.minA[ivar] = blA->data[0];

        report.sumB[ivar] = 0;
        report.maxB[ivar] = blB->data[0];
        report.minB[ivar] = blB->data[0];

        report.sumABDiff[ivar] = 0.;
        report.sumSqVar[ivar] = 0.0;
        
        if(blA->data[0] != blB->data[0])
          report.minError[ivar] = norm(blA->data[0], blB->data[0]);
        else
          report.minError[ivar] = 0;
        /* FIXME should initialize report.maxErrorVal[0][ivar], etc...*/
      }

      cellCount = blA->size[0] * blA->size[1] * blA->size[2];

      /*this is where the acutal grid data is compared*/
      failed = 0;
      
      for(i=0; i<cellCount; i++) {
        if (blA->data[i] != blB->data[i]) {
          if(!failed) {
            failed = 1;
            report.badBlocks[ivar]++;
            if(opts.verbose)
              printf(" %d", blockA);
          }
          

          error = norm(blA->data[i], blB->data[i]);
          /* get sum of square of difference for L2 norm */
          report.sumABDiff[ivar] += pow((fabs(blA->data[i] - blB->data[i])), (double)opts.norm_order);
          if (opts.benchmark=='A') {
            report.sumSqVar[ivar] += pow(fabs(blA->data[i]), (double)opts.norm_order);
          } else {
            report.sumSqVar[ivar] += pow(fabs(blB->data[i]), (double)opts.norm_order);
          }

          if(error > report.maxError[ivar]){
            report.maxError[ivar] = error;
            report.maxErrorVal[0][ivar] = blA->data[i];
            report.maxErrorVal[1][ivar] = blB->data[i];
          }

          if(error < report.minError[ivar]){
            report.minError[ivar] = error;
            report.minErrorVal[0][ivar] = blA->data[i];
            report.minErrorVal[1][ivar] = blB->data[i];
          }

          error = fabs(blA->data[i] - blB->data[i]);
          if(error > report.absError[ivar]){
            report.absError[ivar] = error;
            report.absErrorVal[0][ivar] = blA->data[i];
            report.absErrorVal[1][ivar] = blB->data[i];
          }
        } else { /* we have (blA->data[i] == blB->data[i]) */
          if (opts.benchmark=='A') {
            report.sumSqVar[ivar] += pow(fabs(blA->data[i]), (double)opts.norm_order);
          } else {
            report.sumSqVar[ivar] += pow(fabs(blB->data[i]), (double)opts.norm_order);
          }

        }

        report.sumA[ivar] += blA->data[i]; /* FIXME (volume-weight?) */
        report.maxA[ivar] = (blA->data[i] > report.maxA[ivar]) ? blA->data[i] :
                                                             report.maxA[ivar];
        report.minA[ivar] = (blA->data[i] < report.minA[ivar]) ? blA->data[i] :
                                                             report.minA[ivar];

        report.sumB[ivar] += blB->data[i]; /* FIXME (volume-weight?) */
        report.maxB[ivar] = (blB->data[i] > report.maxB[ivar]) ? blB->data[i] :
                                                             report.maxB[ivar];
        report.minB[ivar] = (blB->data[i] < report.minB[ivar]) ? blB->data[i] :
                                                             report.minB[ivar];

      }

      FR_DeleteBlock(blA);
      FR_DeleteBlock(blB);

    }    

    /* calculate magError */
    error = fabs(report.maxA[ivar]);
    error = fabs(report.minA[ivar]) > error ? fabs(report.minA[ivar]) : error;
    error = fabs(report.maxB[ivar]) > error ? fabs(report.maxB[ivar]) : error;
    error = fabs(report.minB[ivar]) > error ? fabs(report.minB[ivar]) : error;
    error = error > 1.e-99 ? error : 1.e-99;
    report.magError[ivar] = report.absError[ivar] / error;

    /*calc Ln for variable*/
    report.L2Norm[ivar] = pow(report.sumABDiff[ivar], 1.0/(double)opts.norm_order);
    r = pow(report.sumSqVar[ivar], 1.0/(double)opts.norm_order);
    report.L2Norm[ivar] = report.L2Norm[ivar] / 
      ((r > 1.0e-99) ? r : 1.0e-99);

    if (report.magError[ivar] > opts.err_tol) {
      report.fail = 1;
    }

  }

  if(opts.verbose)
    printf("\n");

  /* Create report */
  printf("\n");
  if(report.extraBlocksA)
    printf("%s has %d leaf blocks that don't exist in %s\n", A->filename,
           report.extraBlocksA, B->filename);
  if(report.extraBlocksB)
    printf("%s has %d leaf blocks that don't exist in %s\n", B->filename,
           report.extraBlocksB, A->filename);
  printf("Total leaf blocks compared: %d (all other blocks are ignored)\n", report.blocksCompared);


  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  printf("Var  | Bad Blocks | Min Error ||             Max Error             ||             Abs Error             |\n");
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  printf("     |            |           ||   Error   |     A     |     B     ||   Error   |     A     |     B     |\n");
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  for(i=0; i<report.nvar; i++) {
    printf("%-4s | %-10d | %-9.4g || %-9.4g | %- 9.3g | %- 9.3g || %-9.4g | %- 9.3g | %- 9.3g |\n", 
           report.varnames[i], report.badBlocks[i],

           report.minError[i],

           report.maxError[i],
           report.maxErrorVal[0][i], report.maxErrorVal[1][i],

           report.absError[i],
           report.absErrorVal[0][i], report.absErrorVal[1][i]
           );
  }
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");

  printf("\n");

  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  printf("Var  | L%c-ErrNorm | Mag Error ||                  A                ||                  B                |\n",
   normChar);
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  printf("     | Rel / %c    |           ||    Sum    |    Max    |    Min    ||    Sum    |    Max    |    Min    |\n",
  opts.benchmark);
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
  for(i=0; i<report.nvar; i++) {
    printf("%-4s | %-10.4g | %-9.4g || %- 9.3g | %- 9.3g | %- 9.3g || %- 9.3g | %- 9.3g | %- 9.3g |\n", 
           report.varnames[i], report.L2Norm[i],
           report.magError[i],
           report.sumA[i], report.maxA[i], report.minA[i],
           report.sumB[i], report.maxB[i], report.minB[i]
           );
  }
  printf("-----+------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");


  if (particleReport.particlesCompared > 0) {
    printf("\n");
    printf("\nTotal particles compared: %d\n", particleReport.particlesCompared);

    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    printf("Real Property             | Bad Particles  | Min Error ||             Max Error             ||             Abs Error             |\n");
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    printf("                          |                |           ||   Error   |     A     |     B     ||   Error   |     A     |     B     |\n");
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    for(i=0; i<particleReport.numRealPartProps; i++) {
      printf("%-25s | %-14d | %-9.4g || %-9.4g | %- 9.3g | %- 9.3g || %-9.4g | %- 9.3g | %- 9.3g |\n", 
             particleReport.realPartPropNames[i], particleReport.badParticles[i],

             particleReport.minError[i],

             particleReport.maxError[i],
             particleReport.maxErrorVal[0][i], particleReport.maxErrorVal[1][i],

             particleReport.absError[i],
             particleReport.absErrorVal[0][i], particleReport.absErrorVal[1][i]
             );
    }
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");

    printf("\n");

    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    printf("Real Property             | Bad Particles  | Mag Error ||                  A                ||                  B                |\n");
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    printf("                          |                |           ||    Sum    |    Max    |    Min    ||    Sum    |    Max    |    Min    |\n");
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");
    for(i=0; i<particleReport.numRealPartProps; i++) {
      printf("%-25s | %-14d | %-9.4g || %- 9.3g | %- 9.3g | %- 9.3g || %- 9.3g | %- 9.3g | %- 9.3g |\n", 
             particleReport.realPartPropNames[i], particleReport.badParticles[i],
             particleReport.magError[i],
             particleReport.sumA[i], particleReport.maxA[i], particleReport.minA[i],
             particleReport.sumB[i], particleReport.maxB[i], particleReport.minB[i]
             );
    }
    printf("--------------------------+----------------+-----------++-----------+-----------+-----------++-----------+-----------+-----------+\n");

  }


  if(report.fail || particleReport.fail)
    printf("FAILURE\n");
  else
    printf("SUCCESS\n");
  
  /* Cleanup */
  FR_close(A);
  FR_close(B);
  free(AtoB);
 
  return 0;
}

