VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F

# Test for machine name
include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<


TARGET = nestbat-2.1
TARGET2 = nestbat-2.2

##################################################################################
## Executable for nestbat - old 2.1 version
OBJECTS = \
	nestbat-2.1.o

$(TARGET): $(OBJECTS)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
##################################################################################
## Executable for nestbat - 2.2 version
OBJECTS2 = \
        mod_xc_global.o \
        mod_za_global.o \
        nestbat-2.2.o

$(TARGET2): $(OBJECTS2)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)
##################################################################################



all: $(TARGET) $(TARGET2)

install : all 
	mkdir -p ../../bin_setup
	cp $(TARGET) $(TARGET2) ../../bin_setup


clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.mod $(TARGET) $(TARGET2)
