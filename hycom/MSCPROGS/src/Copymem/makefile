VPATH = .:Build

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  .c


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/  $(INCLUDE)

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

.F.o:
	cd Build/ ; $(CF90) $(CPPFLAGS) $(FFLAGS) $(F77FLG)  $(INCLUDE) -c  ../$<
.F90.o:
	cd Build/ ; $(CF90) $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE)  -c  ../$<





TARGET=copymem

#################################################################################
OBJECTS = \
		  m_parse_blkdat.o   \
		  mod_rand_io.o   \
		  mod_ice_io.o   \
		  mod_evp_io.o   \
		  p_copymem.o  

$(TARGET): $(OBJECTS) 
	cd Build/ ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
##################################################################################

all : $(TARGET)

install : all
	mkdir -p ../../bin
	cp $(TARGET) ../../bin

clean:
	rm Build/*.o Build/*.mod  $(TARGET)

