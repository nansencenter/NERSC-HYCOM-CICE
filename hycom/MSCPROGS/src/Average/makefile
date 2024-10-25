VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F 


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= $(LIBS) -L../../../lib/ -lhycnersc 

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= $(INCLUDE) -I../../../include/ 
# Uncomment below to dump meanssh field
#CPPFLAGS:= $(CPPFLAGS) -DDUMP_MSSH 

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif


# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c  $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE)  -o $*.o ../$<


TARGET = hycave
TARGET2 = ensave
TARGET3 = hycom_mean
all: $(TARGET) $(TARGET2) $(TARGET3)


###############################################################################
OBJECTS= p_hycave.o 

$(TARGET):  $(OBJECTS)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
###############################################################################
OBJECTS2= p_ensave.o

$(TARGET2):  $(OBJECTS2)
	   cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
		###############################################################################
###############################################################################
OBJECTS3= mod_mean.o bigrid.o blkin.o extrct.o getdat.o putdat.o hycom_mean.o

$(TARGET3):  $(OBJECTS3)
	   cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
		###############################################################################

clean:
	rm ./TMP/*.o TMP/*.mod TMP/*.f TMP/*.f90 $(TARGET) $(TARGET2)

install : all
	mkdir -p ../../bin
	cp $(TARGET) $(TARGET2) $(TARGET3) ../../bin
