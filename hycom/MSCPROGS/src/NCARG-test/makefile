VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f 


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


# Rules for running cpp and updating files in TMP directory
# NB: Should be compiled using real4 since it uses ncar graphics 
# real*4 routines
.F90.o:
	cd ./TMP ; $(NCARGCF90) -c $(CPPFLAGS)  $(FFLAGSR4) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(NCARGCF77) -c $(CPPFLAGS) $(FFLAGSR4) $(F77FLG) $(INCLUDE) -o $*.o ../$<

.f.o:
	cd ./TMP ; $(NCARGCF77) -c $(FFLAGSR4) $(F77FLG) $(INCLUDE) -o $*.o ../$<





TARGET  = cplot
TARGET2 = cvplot

all: $(TARGET)  $(TARGET2)
 

###############################################################################
# To generate argocmp
OBJECTS= \
			colorpalette.o \
		   cplot.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(NCARGLD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
# To generate argocmp
OBJECTS2= colorpalette.o \
		    vvumxy.o \
			 m_rotate.o \
			 cvplot.o 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(NCARGLD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################

install : all
	mkdir -p ../../bin
	cp $(TARGET) $(TARGET2) ../../bin

clean:
	rm ./TMP/*.o TMP/*.mod $(TARGET) $(TARGET2)

