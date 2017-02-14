#if 0
    This file contains the indices of the various attributes used within the Proton
    Imaging unit.
#endif

#if 0
    Proton stuff

    Currently the PROTON_ATTRCOUNT has to match the setting of the RAY_ATTR_COUNT
    and EMPROTON_ATTRCOUNT. This has to be done to circumvent the shortcomings of the
    grid particles unit, which uses sending and destination buffers of fixed size
    determined during inititialization of the grid particle unit (gr_ptInit). If
    PROTON_ATTRCOUNT, EMPROTON_ATTRCOUNT and RAY_ATTR_COUNT do not match, the gr_ptInit
    issues a message and stops the calculations.

    PROTON_ATTRCOUNT : total number of proton attributes

    PROTON_POSX : the x coordinate of the proton
    PROTON_POSY : the y coordinate of the proton
    PROTON_POSZ : the z coordinate of the proton
    PROTON_VELX : the x velocity component of the proton
    PROTON_VELY : the y velocity component of the proton
    PROTON_VELZ : the z velocity component of the proton
    PROTON_TIME : the time the proton spends in the domain during a time step
    PROTON_BLCK : block number identifier of proton
    PROTON_PROC : processor number of proton
    PROTON_TAGS : globally unique proton tag
    PROTON_BEAM : beam number where proton originated from
    PROTON_DETC : target detector number of proton
    PROTON_DGJV : the diagnostic total path magnetic current value of the proton
    PROTON_DGKX : the diagnostic total path Bx value of the proton
    PROTON_DGKY : the diagnostic total path By value of the proton
    PROTON_DGKZ : the diagnostic total path Bz value of the proton

#endif

#define PROTON_ATTRCOUNT  16

#define PROTON_POSX 1
#define PROTON_POSY 2
#define PROTON_POSZ 3
#define PROTON_VELX 4
#define PROTON_VELY 5
#define PROTON_VELZ 6
#define PROTON_TIME 7
#define PROTON_BLCK 8
#define PROTON_PROC 9
#define PROTON_TAGS 10
#define PROTON_BEAM 11
#define PROTON_DETC 12
#define PROTON_DGJV 13
#define PROTON_DGKX 14
#define PROTON_DGKY 15
#define PROTON_DGKZ 16

#if 0
    Screen proton stuff

    SCREEN_ATTRCOUNT : total number of screen proton attributes

    SCREEN_POSX : the screen x coordinate of the screen proton
    SCREEN_POSY : the screen y coordinate of the screen proton
    SCREEN_DETC : detector number of screen proton
    SCREEN_DGJV : the diagnostic total path magnetic current value of the screen proton
    SCREEN_DGKX : the diagnostic total path Bx value of the screen proton
    SCREEN_DGKY : the diagnostic total path By value of the screen proton
    SCREEN_DGKZ : the diagnostic total path Bz value of the screen proton

#endif

#define SCREEN_ATTRCOUNT  7

#define SCREEN_POSX 1
#define SCREEN_POSY 2
#define SCREEN_DETC 3
#define SCREEN_DGJV 4
#define SCREEN_DGKX 5
#define SCREEN_DGKY 6
#define SCREEN_DGKZ 7

#if 0
    Proton beam stuff

    BEAM_STRINGLENGTH  : the maximum character string length to accomodate beam info
    BEAM_GRIDARRAYSIZE : the size of each grid array when retrieving all grid points

#endif

#define BEAM_STRINGLENGTH  20
#define BEAM_GRIDARRAYSIZE 100000

#if 0
    Geometry stuff

    GRID_1DCARTESIAN   : Handle for domain geometry
    GRID_2DCARTESIAN   : Handle for domain geometry
    GRID_3DCARTESIAN   : Handle for domain geometry
    GRID_1DCYLINDRICAL : Handle for domain geometry
    GRID_2DCYLINDRICAL : Handle for domain geometry
    GRID_3DCYLINDRICAL : Handle for domain geometry
    GRID_1DSPHERICAL   : Handle for domain geometry
    GRID_2DSPHERICAL   : Handle for domain geometry
    GRID_3DSPHERICAL   : Handle for domain geometry
    GRID_1DPOLAR       : Handle for domain geometry
    GRID_2DPOLAR       : Handle for domain geometry
    GRID_3DPOLAR       : Handle for domain geometry
#endif

#define GRID_1DCARTESIAN    1
#define GRID_2DCARTESIAN    2
#define GRID_3DCARTESIAN    3
#define GRID_1DCYLINDRICAL  4
#define GRID_2DCYLINDRICAL  5
#define GRID_3DCYLINDRICAL  6
#define GRID_1DSPHERICAL    7
#define GRID_2DSPHERICAL    8
#define GRID_3DSPHERICAL    9
#define GRID_1DPOLAR        10
#define GRID_2DPOLAR        11
#define GRID_3DPOLAR        12
