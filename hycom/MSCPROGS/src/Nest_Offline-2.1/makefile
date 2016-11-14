VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F .H .h  

include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

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


TARGET = nestoffline
TARGET2= checknest
TARGET3= nest_layer_remap
TARGET4= clim_weekly
TARGET5= hint_mercator
TARGET6= vremap_mercator
TARGET7= testremap

##################################################################################
## Executable for creating nest files from daily - restart - snapshot etc files
OBJECTS = \
        m_filenesting.o \
		  m_nearestpoint.o \
		  m_rotate2.o\
		  m_parse_blkdat.o \
		  mod_sigma.o \
		  m_bavg_flds_from_ave.o \
		  p_nest_offline.o 


$(TARGET): $(OBJECTS)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
##################################################################################
## Executable for checking nest files 
OBJECTS2 = m_get_nest_info.o \
			  m_get_nest_record.o\
			  p_check_nest.o 


$(TARGET2): $(OBJECTS2)
	cd TMP ; $(LD) $(LINKFLAGS) $(MATLIBS)  -o ../$(TARGET2) $(OBJECTS2) $(LIBS)
##################################################################################
##################################################################################
## Executable for remapping nest files 
OBJECTS3 = mod_sigma.o \
			  m_layer_remapV1.o \
			  m_layer_remapV4.o\
			  m_parse_blkdat.o  \
			  m_dp0kini.o \
		     m_get_nest_info.o \
			  m_get_nest_record.o \
			  m_layer_mixV1.o \
		     m_read_nest_header.o\
			  m_save_nest_header.o \
           p_nest_layer_remap.o


$(TARGET3): $(OBJECTS3)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3) $(OBJECTS3) $(LIBS)
##################################################################################
## Executable for checking nest files 
OBJECTS4 = m_parse_blkdat.o \
		     p_clim_weekly.o 

$(TARGET4): $(OBJECTS4)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET4) $(OBJECTS4) $(LIBS)
##################################################################################
OBJECTS5 = mod_mercatorgrid.o \
			  m_nearestpoint.o \
			  p_hint_mercator.o

$(TARGET5): $(OBJECTS5)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET5) $(OBJECTS5) $(LIBS) 
##################################################################################
OBJECTS6 = mod_sigma.o \
			  m_parse_blkdat.o \
			  m_layer_remapV2.o  \
		     m_rotate2.o \
		     m_layer_mixV1.o \
			  m_pbavg_from_ssh.o \
			  p_vremap_mercator.o

$(TARGET6): $(OBJECTS6)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET6) $(OBJECTS6) $(LIBS)
##################################################################################
##################################################################################
OBJECTS7 = mod_sigma.o \
			  m_layer_remapV4.o \
			  p_test_remap.o

$(TARGET7): $(OBJECTS7)
	cd TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET7) $(OBJECTS7) $(LIBS)
##################################################################################


all: $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) 


clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.mod $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) 
