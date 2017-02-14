#  Makefile for setup StirFromFile

Simulation += Simulation_data.o

Simulation_initBlock.o : Simulation_data.o \
                         Driver_interface.o Grid_interface.o Eos_interface.o
Simulation_init.o : Simulation_data.o Driver_interface.o RuntimeParameters_interface.o
