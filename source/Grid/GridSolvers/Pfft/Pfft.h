#define PFFT_REAL 360
#define PFFT_REAL2C 365
#define PFFT_REAL2C_STUFF 366
#define PFFT_REAL2C_EXTEND 367
#define PFFT_COMPLEX 768
#define PFFT_COMPLEX_STUFFED 769
#define PFFT_COS 657
#define PFFT_SIN 645
#define PFFT_COSQ 658
#define PFFT_SINQ 646
#define PFFT_COS_CC 659
#define PFFT_SIN_CC 647
#define PFFT_COS_IV 660
#define PFFT_SIN_IV 648
#define PFFT_TRANSFORM_NONE -2

#define PFFT_FORWARD 346
#define PFFT_INVERSE 323

#define PFFT_PCLDATA_REAL 365
#define PFFT_PCLDATA_COMPLEX 768
#define PFFT_PCLDATA_COMPLEX_STUFFED 769
#define PFFT_PCLDATA_COMPLEX_EXTENDED 367

#define PFFT_SPOSI 1
#define PFFT_SPOSJ 2
#define PFFT_SPOSK 3
#define PFFT_EPOSI 4
#define PFFT_EPOSJ 5
#define PFFT_EPOSK 6
#define PFFT_BLKID 7
#define PFFT_BSIZE 8
#define PFFT_MAPEND 8
#define PFFT_NUMBLKS 1

#define PFFT_BND_DIRICHLET 2
#define PFFT_BND_NEUMANN   3
#define PFFT_BND_PERIODIC  1

#define TO_PFFT 1
#define FROM_PFFT 2

#if 0
We provide two definitions which can be used as optional arguments for jProcs and kProcs 
in Grid_pfftInit() call.  These specify that work should be distributed over a single 
dimension only in a 3D PFFT grid.  

These prevent an occasional abort, which may happen if the user provides jProcs and kProcs
but the PFFT framework concludes there is insufficient work to spread among available processors.
This is a very real possibility when PFFT is used with Multigrid!  If you wish to 
distribute over one dimension only, you should use the definitions.  We make the 
values negative to avoid interpreting a users selection wrongly in Grid_pfftInit.
#endif

#define PFFT_ALL_PROCS -100
#define PFFT_ONE_PROC -200


#define PFFT_FLASH_NODE -300
#define PFFT_PENCIL_NODE -400

#define PFFT_SINGLE_LEVEL 123
#define PFFT_MIXED_LEVEL 124
