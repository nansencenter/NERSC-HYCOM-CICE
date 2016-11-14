VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F 


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc $(LIBS) 

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../ -I../../../include/ $(INCLUDE)

# Undefine matlab output (too much hassle from machine to machine, old outdated
# libraries etc etc ..)
CPPFLAGS:= $(CPPFLAGS) -UDUMP_MSSH -UMATLAB

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

TARGET  = ssh_from_restart
all: $(TARGET) 


###############################################################################
# To generate argocmp
OBJECTS=\
	mod_sigma.o\
	mod_read_rstab.o\
	m_parse_blkdat.o \
	p_ssh_from_state.o 

$(TARGET): $(OBJECTS)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################

install : all 
	mkdir -p     ../../bin
	cp $(TARGET) ../../bin

clean:
	rm TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.H TMP/*.mod $(TARGET)
