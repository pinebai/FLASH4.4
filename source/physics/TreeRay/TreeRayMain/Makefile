#	Makefile for TreeRad module
TreeRay += healpix.o healpix_types.o bit_manipulation.o \
           tr_bhLocalInterface.o tr_bhGenIntersectList.o \
           tr_bhRadToGas.o \
           TreeRay_data.o tr_bhFinalizeCell.o


TreeRay += tr_odAccBotNode.o tr_odAccNode.o \
           tr_odBTInit.o tr_odFillBotNode.o tr_odFinalize.o \
           tr_odFinalizeCell.o tr_odGetNodeStruct.o tr_odInit.o \
           tr_odInitFieldVar.o tr_odInterface.o tr_odStartBlock.o \
           tr_odNodeContrib.o tr_odSelfContrib.o 


tr_bhGenIntersectList.o : healpix.o
healpix.o : bit_manipulation.o
bit_manipulation.o : healpix_types.o
