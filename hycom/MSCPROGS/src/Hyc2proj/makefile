VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ -I../ $(INCLUDE)


CPPFLAGS:= $(CPPFLAGS) -UMULT_GRP_PER_FILE 
#CPPFLAGS:= $(CPPFLAGS) 

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90; cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o $*.f90

.F.o:
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./TMP/$*.f; cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o $*.f





TARGET2 = hyc2proj
TARGET3 = hyc2stations

all: $(TARGET2) $(TARGET3)


###############################################################################
# To generate hyc2proj
OBJECTS2=\
	m_nearestpoint.o \
	m_mersea_prepare.o\
	m_process_arguments.o\
	m_fields_to_plot.o\
	mod_toproj.o\
	mod_rotate.o\
	m_strmf_eval.o \
	m_mixlayer_depths.o \
	mod_station.o  \
	mod_netcdf_file.o  \
	mod_levitus.o  \
	m_read_mean_ssh.o \
	m_construct_filename.o \
	m_bio_conversions.o \
	p_hyc2proj.o 

$(TARGET2): $(OBJECTS2)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
# To generate hyc2stations
OBJECTS3= \
	mod_rotate.o  \
	mod_levitus.o \
	mod_station.o  \
	mod_toproj.o\
	mod_netcdf_file.o \
	m_fields_to_plot.o \
	m_construct_filename.o \
	p_hyc2stations.o 

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################

clean:
	rm ./TMP/*.f* ./TMP/*.o TMP/*.h TMP/*.mod $(TARGET2) $(TARGET3) 

install : all
	mkdir -p ../../bin
	cp hyc2proj hyc2stations ../../bin

