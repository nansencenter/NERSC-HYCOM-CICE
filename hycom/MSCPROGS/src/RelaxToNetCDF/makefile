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

TARGET2=hclimlevels
TARGET3=hclimtonc

 
all: $(TARGET2) $(TARGET3)

###############################################################################
# To generate hclimlevels
OBJECTS2 = p_hycomclim_levels.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
# To generate hclimtonc
OBJECTS3 = \
		   m_handle_err.o  \
			m_nearestpoint.o  \
		   mod_toproj.o  \
			p_hclimtonc.o  \
			mapxy.o

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################

install : all 
	mkdir -p ../../bin/
	cp $(TARGET2) $(TARGET3) ../../bin/

clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET2) $(TARGET3)
