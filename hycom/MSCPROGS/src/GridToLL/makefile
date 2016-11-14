VPATH = .:TMP
.SUFFIXES: .F90 .F .o

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

.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(F90FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(F77FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<


TARGET  = gridtoll
TARGET2  = findnearest

OBJECTS = m_handle_err.o\
			 p_gridtoll.o 
OBJECTS2 = m_handle_err.o\
			 m_spherdist.o\
			 p_findnearest.o 



all: $(TARGET)  $(TARGET2)

$(TARGET):  $(OBJECTS) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)

$(TARGET2):  $(OBJECTS2) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)


clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET) $(TARGET2)



