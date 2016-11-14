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



# Set wether to use TRIP or TRIP05 data base and ERA-I or ERA40
CPPFLAGS:=$(CPPFLAGS) -DTRIP05 -UTRIP 


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


TARGET = trip_weights
TARGET2 = trip_paths
TARGET3 = trip_flow
TARGET4 = trip_tohycom
TARGET5 = trip_tobiorivGN

 
all: $(TARGET)  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5)


###############################################################################
OBJECTS=m_handle_err.o \
		  mod_trip.o \
	m_read_runoff_era40.o \
	m_read_runoff_erai.o \
	p_trip_riverweights.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
OBJECTS2=\
	mod_trip.o \
	p_trip_riverpaths.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
OBJECTS3=m_handle_err.o \
		  mod_trip.o \
			m_read_runoff_era40.o \
			m_read_runoff_erai.o \
			m_era40_fix.o \
			p_trip_riverflow.o

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################
OBJECTS4=m_handle_err.o \
			m_nearestpoint.o \
			mod_trip.o \
			p_trip_tohycom.o

$(TARGET4): $(OBJECTS4) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET4)  $(OBJECTS4) $(LIBS) 
###############################################################################
OBJECTS5=m_handle_err.o \
			m_nearestpoint.o \
			mod_trip.o \
			p_trip_tobiorivGN.o

$(TARGET5): $(OBJECTS5) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET5)  $(OBJECTS5) $(LIBS) 
###############################################################################


clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET)  $(TARGET2)  $(TARGET3)  $(TARGET4) $(TARGET5)

install: all
	mkdir -p ../../bin_setup/
	cp $(TARGET) $(TARGET2)  $(TARGET3)  $(TARGET4)  ../../bin_setup/

