VPATH = .:Build/

.SUFFIXES:
.SUFFIXES: .o .F90 .F 


# Tests dont work with unix make ? (works with gmake ...)
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

# Rules for running cpp and updating files 
.F.o:
	cd Build; $(CF90) $(CPPFLAGS) $(INCLUDE) -c $(FFLAGS) $(F77FLG) -o $@ ../$<

.F90.o:
	cd Build; $(CF90) $(CPPFLAGS) $(INCLUDE) -c $(FFLAGS) $(F90FLG) -o $@ ../$<


TARGET=forfun_nersc-2.1

all: $(TARGET) 


#################################################################################
OBJECTS = \
		  m_era40_fix.o \
		  m_bilin_ecmwf2.o  \
		  m_bilin_ncep_gauss.o \
		  m_ncvar_dims.o  \
		  m_ncvar_read.o \
		  mod_forcing_nersc.o  \
		  m_rotate.o \
		  m_rotate_sphere.o \
		  m_relhumid.o \
		  m_satvap.o \
		  m_qsw0.o \
		  m_read_ecmwf_nc.o \
		  m_read_ecmwf.o  \
		  m_read_ncep.o \
		  m_read_clim.o \
		  m_read_era40.o \
		  mod_storm.o \
		  p_forfun_nersc_new.o 

$(TARGET): $(OBJECTS) 
	cd Build; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
#################################################################################


install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) ../../bin_setup

clean:
	rm Build/*.o Build/*.mod  $(TARGET)

