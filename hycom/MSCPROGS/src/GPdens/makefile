VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F   

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
	cd ./TMP ; $(CF90) $(CPPFLAGS) -c $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) $(CPPFLAGS) -c $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<



TARGET1 = gpdens
all: $(TARGET1) 

#############################################################################
OBJECTS1= \
			 mod_data_new.o \
			 mod_read_gp.o \
			 mod_netcdf_ops.o \
			 m_exceedence.o \
			 m_mkestat.o \
			 m_correlation.o \
			 m_accplot.o \
			 m_mdian1.o \
			 m_shfact.o \
			 m_shfilt.o \
			 m_time_filter.o \
			 m_density.o \
			 m_densitylevels.o \
			 m_rotate.o \
			 p_gpdens.o

$(TARGET1): $(OBJECTS1) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET1) $(OBJECTS1) $(LIBS) 
#############################################################################

install: all
	mkdir -p ../../bin
	cp $(TARGET1) ../../bin


clean:
	rm ./TMP/*.o  TMP/*.mod $(TARGET1)
