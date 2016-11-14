VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F  



#CPP id no longer used
include ../Make.Inc/make.inc


# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc  $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/  -I../ $(INCLUDE)

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif


# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c  $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS)  $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<


TARGET = igrid
TARGET2= irstrt


all: $(TARGET) $(TARGET2)

############################################################################
OBJECTS=\
		  m_parse_blkdat.o\
		  m_griddump.o\
		  igrid.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
############################################################################
OBJECTS2=\
			m_parse_blkdat.o\
			mod_sigma.o\
			irstrt.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
############################################################################

clean:
	rm $(TARGET2) $(TARGET)
	cd ./TMP ; rm *.o *.f90 *.h *.mod

