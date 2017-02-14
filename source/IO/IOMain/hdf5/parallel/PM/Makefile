IO += io_writeData.o io_readData.o 


io_writeData.o : RuntimeParameters_getAll.o IO_getScalar.o Logfile_stamp.o
io_readData.o :RuntimeParameters_setPrev.o IO_getScalar.o
IO_writeCheckpoint.o : io_writeData.o 


