#!/bin/bash -x
set -e

HDF5_PATH='/gpfs/vesta_scratch/projects/TurbNuclComb_esp/cdaley/vesta/software/V1R2M0/hdf5-fen/1.8.10/xl/ibm-compilers-nov2012'
PATCH_FILE='no_particles_comparison.diff'
OLD_DIR=$(pwd)

cd ../../tools/sfocu
make clean

# Create a regular build of sfocu.  Rename the binary to make it clear
# that it will compare particles.  The particle comparison is
# extremely slow (N^2 algorithm)!
echo -e "\nDo a regular sfocu build..."
make NO_MPI=True NO_NCDF=True CCOMP=xlc_r CFLAGS_HDF5="-I${HDF5_PATH}/include -DH5_USE_16_API" LIB_HDF5="-L${HDF5_PATH}/lib -lhdf5"

mv sfocu sfocu.particles

echo -e "\nDo a custom sfocu build (without the particle comparison)"
# Create a custom build of sfocu.  Remove the particle comparison.
cat > ${PATCH_FILE} << EOF
Index: flash_reader_hdf5.c
===================================================================
--- flash_reader_hdf5.c	(revision 18983)
+++ flash_reader_hdf5.c	(working copy)
@@ -509,7 +509,7 @@
       out->varnames[i][FR_VAR_STRING_SIZE] = '\0';
     }
   }
-  FR_GetNumParticles_HDF5(&handle, &num_particles);
+  FR_GetNumParticles_HDF5(&handle, &num_particles); num_particles = 0;
 
   out->totalparticles = num_particles;
   out->numRealPartProps = 0;
EOF

patch -p0 < ${PATCH_FILE}
make NO_MPI=True NO_NCDF=True CCOMP=xlc_r CFLAGS_HDF5="-I${HDF5_PATH}/include -DH5_USE_16_API" LIB_HDF5="-L${HDF5_PATH}/lib -lhdf5"

# Clean up
patch -R -p0 < ${PATCH_FILE}
rm ${PATCH_FILE}
rm -f *.o
cd ${OLD_DIR}
