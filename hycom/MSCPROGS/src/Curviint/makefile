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


TARGET = curviint-2.1
TARGET2= curviint-2.2
TARGET3= curviint_sed-2.2

 
all: $(TARGET)  $(TARGET2)  $(TARGET3)

###############################################################################
OBJECTS=\
	mod_xc_global.o\
	mod_xc_local.o\
	mod_za_global.o\
	mod_za_local.o \
	p_curviint-2.1.o 

$(TARGET): $(OBJECTS)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
OBJECTS2=\
	m_parse_blkdat.o\
	mod_xc_global.o\
	mod_xc_local.o\
	mod_za_global.o\
	mod_za_local.o \
	mod_curviint_interp.o \
	p_curviint-2.2.o 

$(TARGET2): $(OBJECTS2)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
OBJECTS3=\
	m_parse_blkdat.o\
	mod_xc_global.o\
	mod_xc_local.o\
	mod_za_global.o\
	mod_za_local.o \
	mod_curviint_interp.o \
	p_curviint_sed-2.2.o 

$(TARGET3): $(OBJECTS3)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################
install : all
	mkdir -p     ../../bin_setup
	cp $(TARGET) $(TARGET2) $(TARGET3) ../../bin_setup

clean:
	rm ./TMP/*.o TMP/*.mod $(TARGET) $(TARGET2)
