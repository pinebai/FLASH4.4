#	Makefile for gravity module:  plane-parallel gravitational field

#	Additional dependencies

Gravity += Gravity_data.o grv_accelOneZone.o

Gravity_accelOneRow.o  : grv_accelOneZone.o Gravity_data.o
grv_accelOneZone.o : Gravity_data.o
Gravity_init.o         : Gravity_data.o
