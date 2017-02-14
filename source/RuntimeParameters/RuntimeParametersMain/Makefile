RuntimeParameters += RuntimeParameters_data.o \
	rp_getOpt.o \
	rp_getArgument.o \
	rp_initParameters.o \
	rp_rules.o rp_storeIgnoredParams.o


RuntimeParameters_getAll.o : RuntimeParameters_data.o nameValueLL_getAll.o 
RuntimeParameters_data.o : nameValueLL_data.o
RuntimeParameters_init.o : RuntimeParameters_data.o
RuntimeParameters_get.o : RuntimeParameters_data.o nameValueLL_get.o 
RuntimeParameters_getPrev.o : RuntimeParameters_data.o nameValueLL_get.o 
RuntimeParameters_setPrev.o : RuntimeParameters_data.o nameValueLL_set.o 
RuntimeParameters_set.o : RuntimeParameters_data.o nameValueLL_set.o 
RuntimeParameters_read.o : RuntimeParameters_set.o
RuntimeParameters_getNum.o : RuntimeParameters_data.o nameValueLL_getNum.o
rp_initParameters.o : RuntimeParameters_add.o 
RuntimeParameters_stampIgnored.o : Logfile_stampMessage.o

.SECONDARY: f90_unix_env.mod
