VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F 

include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS) 

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)

CPPFLAGS:=$(CPPFLAGS) -DNORSEXGRID


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


TARGET  = norsexclim

 
all: $(TARGET) 


###############################################################################
# To generate hyc2proj
OBJECTS=mod_toproj.o\
		  m_handle_err.o\
		  p_norsexclim.o\
		  norsex_locate.o\
		  norsex_mapxy.o\
		  norsex_mapll.o \
        norsex_getgrid.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################

install : all
	mkdir -p     ../../bin
	cp $(TARGET) ../../bin

clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET)

