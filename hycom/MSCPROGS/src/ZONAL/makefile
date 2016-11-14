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



TARGET  = zonal
TARGET2 = mosf

all : $(TARGET) 

########################################################################
OBJECTS= \
	m_read_blkdat.o \
	m_fields_to_plot.o \
	m_handle_err.o \
	p_zonal.o 

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#############################################################################

OBJECTS2= \
	mod_overturning.o \
	m_read_blkdat.o \
	m_fields_to_plot.o \
	m_handle_err.o \
	p_mosf.o 


$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
#############################################################################

all : $(TARGET) $(TARGET2)

# We dont trust mosf yet
install : all 
	mkdir -p ../../bin
	cp $(TARGET) ../../bin


clean:
	rm TMP/*.o TMP/*.mod $(TARGET) $(TARGET2)
