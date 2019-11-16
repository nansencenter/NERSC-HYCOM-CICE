VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F  


include ../../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../../lib -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../../include/ $(INCLUDE)


# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif


# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c  $(FFLAGS) $(F90FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c  $(FFLAGS) $(F77FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<


TARGET = grid


all: $(TARGET)

###########################################################################
OBJECTS = \
        mod_sphere_tools.o\
        m_spherdist.o\
        m_avecell.o\
        m_bigrid.o\
        m_pivotsearch.o\
        m_bilin.o\
        m_bilinibcao.o\
        m_checktopo.o\
        m_grid_consistency.o\
        m_quick_sort.o\
        m_ibcao_median.o\
        m_invdistibcao3.o\
        m_invdistibcao4.o\
        mod_bodc.o\
        mod_conman.o\
        mod_etopo5.o\
        mod_gebco.o\
        mod_ibcao.o\
        m_sfilt.o\
        m_shapiro_filt.o\
        m_shapiro_ini.o\
        m_tecconfgrid.o\
        m_tecfld.o\
        m_ncfld.o\
        m_topofix.o\
        m_topography.o \
        m_grid_to_hycom.o \
        felt_inout.o\
        grid.o 

#        mod_confmap.o\
#        m_oldtonew.o\
#        m_pivotp.o\
#        m_newtoold.o\
#        m_initconfmap.o\

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
###########################################################################

all : $(TARGET)

install : all
	mkdir -p ../../../bin_setup/
	cp $(TARGET) ../../../bin_setup/

clean:
	rm TMP/*.o TMP/*.mod TMP/*.MOD $(TARGET)

