
!!#define DEBUG_SOLID 1
#define WRITEFORCES 1
#define WRITESTATES 1

#define TRIANG_MAPPING_A

!! Element Definitions from Gmsh Meshing software
!! standard. Look in http://geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
!! If you have another input ordering translate to this!
#define TWO_NODE_LINE            1
#define THREE_NODE_TRIANGLE      2
#define FOUR_NODE_QUADRILATERAL  3
#define FOUR_NODE_TETRAHEDRON    4
#define EIGHT_NODE_HEXAHEDRON    5
#define THREE_NODE_LINE          8
#define SIX_NODE_TRIANGLE        9
#define NINE_NODE_QUADRILATERAL 10
#define TEN_NODE_TETRAHEDRON    11
#define TWSEVEN_NODE_HEXAHEDRON 12  
#define EIGHTEEN_NODE_PRISM     13
#define FOURTEEN_NODE_PYRAMID   14
#define ONE_NODE_POINT          15


!! Node Numbers
#define ONE_NODE    1
#define TWO_NODES   2
#define THREE_NODES 3
#define FOUR_NODES  4
#define FIVE_NODES  5
#define SIX_NODES   6
#define EIGHT_NODES 8
#define NINE_NODES  9
#define TEN_NODES  10
#define THIRTEEN_NODES 13
#define FOURTEEN_NODES 14
#define FIFTEEN_NODES  15
#define EIGHTEEN_NODES 18
#define TWENTY_NODES   20
#define TWENTYSEVEN_NODES 27

!! BodyType Definitions
#define BODYTYPE_RIGID      1
#define BODYTYPE_2DFLEXIBLE 2
#define BODYTYPE_3DFLEXIBLE 3
#define BODYTYPE_RBC        4

!! Restraint Definitions
#define NDMAX 3
#define NMAXPARAMRES 10

!! Material Behavior
#define MATERIAL_KIRCHHOFF      1
#define MATERIAL_BIOT           2
#define MATERIAL_WLC_POW        3
#define MATERIAL_WLC_CQ         4
#define MATERIAL_WLC_FENE       5
#define MATERIAL_SF_WLC_POW     6
#define MATERIAL_SF_WLC_CQ	7
#define MATERIAL_SF_WLC_FENE    8
#define MATERIAL_AREA_VOL	9
#define MATERIAL_SF_AREA_VOL    10
#define MATERIAL_BEND           11
#define MATERIAL_SF_BEND        12
#define MATERIAL_FENE_POW       13
#define MATERIAL_SF_FENE_POW    14

!! Time-marching Integrator
#define SOLIDINTEG_GENALPHA     1
#define SOLIDINTEG_PREDCORR     2
#define SOLIDINTEG_MODVVERLET   3

!! Matrix Storage Scheme (only CSC is tested)
#define FEM_CSR 1
#define FEM_CSC 2
#define	FEM_MATFORMAT 2

! rbc_visc models
#define Non_viscous    0
#define visc_relVel1   1
#define visc_relVel2   2
#define Pep_viscModel  3

! Define SolidMechanics Procedure flags
#define SM_INIT          11
#define SM_ADVANCE       12
#define SM_ADVANCE1DT    13
#define SM_CHECKCONVERG  14
#define SM_FINALIZE      15
#define SM_WRITECHECKPT  16
#define SM_TRUE          1
#define SM_FALSE         0

! Define SM iopt flags
! These can be defined for an individual body type (so duplicates are ok)
#define SM_IOPT_QN 1
#define SM_IOPT_QI 2
#define SM_IOPT_NSTEP 3
#define SM_IOPT_NMIDSTEP 4

! Initial Conditions Flags
! 0= just assume that the ICs are all zero
#define SM_IC_ZERO     0 
! 1= load an hdf5 file per body that defines the ICs
#define SM_IC_EXTFILE  1 
! 2= apply rigid body motion from Prescribed Kinematics
#define SM_IC_PK       2 

! restraint flags
#define ALL_DOF 0

!! Corrections converged
#define SM_NOTCONVERGED 0
#define SM_CONVERGED    1
#define SM_TESTCONVERGE    121
#define SM_SETNOTCONVERGED 122
 
#define MAXNODERBC 4

!! Writing options to Tecplot
#define WRITEPOS    1
#define WRITEVEL    2
#define WRITEFORC   3
#define WRITEHDF5   4

!! Rigid Body Transformations
#define BODYFRAME_NODE 1
#define RB_IDENTITY 1000
#define RB_EULER321 1321
#define RB_QUATERNN 1111 
#define RB_TWODIM   1002

!! Rigid Body Analityical types
#define RB_MAXANNPARAM   4   
#define RB_ANNELLIPSOID 54
#define RB_ANNSPHERE    55
#define RB_ANNDISC      56
#define RB_ANNRBC       57

!! Angle location for Euler angle sequences
#define ANGLE_1  1
#define ANGLE_2  2
#define ANGLE_3  3
