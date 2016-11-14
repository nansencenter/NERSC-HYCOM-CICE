VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  

include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc $(LIBS)

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


TARGET = ensstat_field
TARGET2= ensstat_point
scripts= ensstat_field3D.sh  ensstat_point3D.sh
all : $(TARGET) $(TARGET2)

#############################################################################
OBJECTS= \
			m_parse_blkdat.o \
         p_ensstat_field.o 

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#############################################################################
OBJECTS2= \
			m_parse_blkdat.o \
         p_ensstat_point.o 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
#############################################################################

install : all
	mkdir -p ../../bin
	cp $(TARGET) $(TARGET2) $(scripts) ../../bin

clean:
	rm TMP/*.o TMP/*.mod $(TARGET) $(TARGET2)
