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

CPPFLAGS:= $(CPPFLAGS) -UMATLAB

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c  $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c  $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<


TARGET2= icedrift2
TARGET3= icedrift_osisaf
targets = $(TARGET2) $(TARGET3)
script=gen_icedrift_ens.sh
all : $(TARGET2) $(TARGET3)


#############################################################################
OBJECTS2= \
	m_ncvar_dims.o\
	m_ncvar_read.o\
	m_rk2.o \
	p_icedrift2.o 

OBJECTS3= \
	m_ncvar_dims.o\
	m_ncvar_read.o\
	m_rk2.o\
	p_icedrift_osisaf.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS) 
$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3) $(OBJECTS3) $(LIBS) 
#############################################################################
		 
install : all
	mkdir -p ../../bin
	cp $(TARGET2) $(TARGET3) $(script) ../../bin

clean:
	rm ./TMP/*.o TMP/*.mod $(TARGET2) $(TARGET3)


