VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIB_NC)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INC_NC)



# Set wether to use TRIP or TRIP05 data base
#CPPFLAGS:=($CPPFLAGS) -UTRIP05 -DTRIP
CPPFLAGS:=$(CPPFLAGS) -DTRIP05 -UTRIP


# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<


TARGET = trip_weights
TARGET2 = trip_paths
TARGET3 = trip_riverflow
TARGET4 = trip_tohycom

 
all: $(TARGET)  $(TARGET2) $(TARGET3)


###############################################################################
OBJECTS=m_handle_err.o \
	mod_trip.o \
	m_read_runoff_era40.o \
	p_trip_riverweights.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
OBJECTS2=p_trip_riverpaths.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
OBJECTS3=m_handle_err.o \
			m_read_runoff_era40.o \
			p_trip_riverflow.o

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################
OBJECTS4=m_handle_err.o \
			p_trip_tohycom.o

$(TARGET4): $(OBJECTS4) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET4)  $(OBJECTS4) $(LIBS) 
###############################################################################


clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod

