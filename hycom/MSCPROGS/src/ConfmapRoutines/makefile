VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F  


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib -lhycnersc -lconfmap $(LIBS)

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
	cd ./TMP ; $(CF90) -c  $(FFLAGS) $(F90FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c  $(FFLAGS) $(F77FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<


TARGET1 = confmap_ll2ij
TARGET2 = confmap_gmap
all : $(TARGET1) $(TARGET2)

###########################################################################
OBJECTS1= \
			 p_ll2ij.o

OBJECTS2= \
			 p_confmap_gmap.o

$(TARGET1): $(OBJECTS1) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET1) $(OBJECTS1) $(LIBS) 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
###########################################################################

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) ../../bin_setup/

clean:
	rm ./TMP/*.o TMP/*.mod TMP/*.MOD $(TARGET)

