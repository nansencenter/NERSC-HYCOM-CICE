VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F .H   

include ../Make.Inc/make.inc


# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../ -I../../../include/ $(INCLUDE)


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



TARGET  = thkdf4-2.2.37
TARGET2  = sssrmx-2.2.37

all : $(TARGET) $(TARGET2)

########################################################################
OBJECTS= \
	p_thkdf4.o 

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#############################################################################
OBJECTS2= \
	p_sssrmx.o 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
#############################################################################


# We dont trust mosf yet
install : all 
	mkdir -p ../../bin
	cp $(TARGET) ../../bin


clean:
	rm TMP/*.o TMP/*.mod $(TARGET) $(TARGET2)
