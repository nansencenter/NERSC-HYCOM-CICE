VPATH = .:TMP
.SUFFIXES:
.SUFFIXES: .F90 .F .o .c

include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ 


# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(F90FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(F77FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<



TARGET = vint 
TARGET2= hint

OBJECTS = p_vint.o 
OBJECTS2= p_hint.o 



all: $(TARGET)  $(TARGET2)

$(TARGET):  $(OBJECTS) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)

$(TARGET2):  $(OBJECTS2) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)

#install: all
#	mkdir -p ../../bin/
#	cp $(TARGET) ../../bin


clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod



