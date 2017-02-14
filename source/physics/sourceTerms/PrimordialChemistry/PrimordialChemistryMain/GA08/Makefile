#	Makefile for PrimordialChemistry module

PrimordialChemistry += PrimordialChemistry_data.o  PrimordialChemistry_dataEOS.o \
	     pchem_interface.o pchem_burner.o pchem_azbar.o \
             pchem_mapNetworkToSpecies.o pchem_initNetwork.o \
             pchem_dataNetworkSize.o pchem_network.o pchem_networkRates.o \
 	     pchem_networkScreen.o pchem_networkSparsePointers.o \
   	     pchem_networkDenseJakob.o pchem_gift.o pchem_data.o  \
             pchem_matrix.o pchem_coolFunction.o pchem_coolMole2.o

PrimordialChemistry_data.o : pchem_dataNetworkSize.o
pchem_networkSparseJakob.o : PrimordialChemistry_data.o pchem_data.o