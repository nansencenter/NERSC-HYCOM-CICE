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
	cd ./TMP ; $(CF90) -c $(FFLAGS) $(CPPFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(FFLAGS) $(CPPFLAGS)  $(F77FLG) $(INCLUDE) -o $*.o ../$<



TARGET =  csr2mod_GE

#################################################################################
# KAL - put these in the right order (depends on use module), and dependencies shouldnt be necessary
OBJECTS = m_admit2.o \
          m_constituents.o \
          m_lpeqmt.o \
          m_tptid1.o \
          m_tptide.o \
          m_csrtide.o \
          m_har_mdomain.o \
          m_harmonics.o \
          m_julday.o \
          m_kalday.o \
			 csr2mod.o 


$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#################################################################################


all : $(TARGET)

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) ../../bin_setup/

clean:
	rm  TMP/*.o TMP/*.mod $(TARGET)

