
#ifdef REORDERED
#define PM_REORDERED
#else
#undef PM_REORDERED
#endif

#ifdef PM_REORDERED

#define ARRAY_FIVE(arr,i,x,y,z,b) arr(x,y,z,i,b)
#define ARRAY_FOUR(arr,i,x,y,z) arr(x,y,z,i)
#define ARRAY_THREE(arr,x,y,z) arr(x,y,z)

#else

#define ARRAY_FIVE(arr,i,x,y,z,b) arr(i,x,y,z,b)
#define ARRAY_FOUR(arr,i,x,y,z) arr(i,x,y,z)
#define ARRAY_THREE(arr,x,y,z) arr(x,y,z)

#endif

