# Makefile for source/physics/sourceTerms/Burn/nuclearBurn/Aprox13

Burn += bn_mapNetworkToSpecies.o \
	bn_dataAprox13.o bn_initNetwork.o bn_dataNetworkSize.o \
	bn_network.o bn_networkRates.o bn_networkScreen.o bn_networkTable.o \
	bn_networkSparsePointers.o bn_networkDenseJakob.o bn_networkSparseJakob.o \
	bn_gift.o 

Burn_data.o : bn_dataNetworkSize.o
bn_networkSparseJakob.o : Burn_data.o bn_dataAprox13.o 

