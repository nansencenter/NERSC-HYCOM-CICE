VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F  .f


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc  $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../ -I../../../include/  $(INCLUDE)


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

.f.o:
	cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) $(INCLUDE)  -o $*.o ../$<

TARGET = rmunew
TARGET2 = z_woa2013_bio

all: $(TARGET) $(TARGET2)


###############################################################################
OBJECTS= rmunew.o \
			zh.o

$(TARGET):  $(OBJECTS)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
OBJECTS2= z_woa2013_bio.o \
			zh.o \
			interp.o

$(TARGET2):  $(OBJECTS2)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################

clean:
	rm ./TMP/*.o TMP/*.mod  $(TARGET) $(TARGET2) 

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) ../../bin_setup/
	cp $(TARGET2) ../../bin_setup/
