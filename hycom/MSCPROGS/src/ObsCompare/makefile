VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

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
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(CPPFLAGS) $(INCLUDE) -o $*.o ../$<

TARGET  = argocmp
TARGET2 = argoreduce
TARGET3 = argo_mersea
TARGET4 = slacmp
TARGET5 = mldcmp
scripts= argo_extract_profiles_new.sh argo_extract_profiles.sh slacmp_ftp.sh

 
all: $(TARGET)  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5)

###############################################################################
# To generate argocmp
OBJECTS=  \
			 mod_levitus.o \
			 m_handle_err.o \
          m_read_coriolis.o \
			 m_dump_coriolis.o \
		    p_argocmp.o  

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
# To generate argoreduce
OBJECTS2= \
			 mod_sphere_tools.o \
			 mod_regions.o \
			 p_argoreduce.o 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
# To generate argo_mersea
OBJECTS3= p_argo_mersea.o 

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################
# To generate slacmp
OBJECTS4=\
			 mod_sphere_tools.o \
			 mod_regions.o \
          m_duacsread.o \
			 p_slacompare.o 

$(TARGET4): $(OBJECTS4) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET4)  $(OBJECTS4)  $(LIBS)
###############################################################################
###############################################################################
OBJECTS5= \
          mod_netcdf_helpers.o \
			 mod_sphere_tools.o \
          mod_regions.o \
          m_dBM_MLD_read.o \
			 p_mldcompare.o 

$(TARGET5): $(OBJECTS5) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET5)  $(OBJECTS5)  $(LIBS)
###############################################################################


install : all
	mkdir -p ../../bin
	cp $(TARGET)  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) ../../bin
	cp $(scripts) ../../bin


clean:
	rm TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET)  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) 



