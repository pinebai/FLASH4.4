../setup -auto PreHeatPlug -2d +cylindrical -nxb=16 -nyb=16 +hdf5typeio \
species=cham,targ,gas,LEH,wash +mtmmmt +laser +uhd3t \
+mgd mgd_meshgroups=12 \
-parfile=flash.par -maxblocks=500 \
-objdir=ImplosionPlug_1 -makefile=gnu -noclobber
