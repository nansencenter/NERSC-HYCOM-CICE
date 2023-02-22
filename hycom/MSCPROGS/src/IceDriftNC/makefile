VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F 

include ../Make.Inc/make.inc


# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)

#CPPFLAGS:= $(CPPFLAGS) -UMATLAB

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90; cd ./TMP ; $(CF90) -c  $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./TMP/$*.f; cd ./TMP ; $(CF77) -c  $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<


TARGET0= icedriftnc_osisaf
targets = $(TARGET0)
script=gen_icedrift_ens.sh
all : $(TARGET0)


#############################################################################
OBJECTS0= \
	m_ncvar_dims.o\
	m_ncvar_read.o\
	m_rk2.o\
	m_handle_err.o\
	m_ncvar_dump.o\
	p_icedrift_osisaf.o

$(TARGET0): $(OBJECTS0) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET0) $(OBJECTS0) $(LIBS) 
#############################################################################
		 
install : all
	mkdir -p ../../bin
	cp $(TARGET0) $(script) ../../bin

clean:
	rm ./TMP/*.o TMP/*.mod $(TARGET0)


