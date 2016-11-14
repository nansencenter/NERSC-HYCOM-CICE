VPATH = .:TMP
.SUFFIXES:
.SUFFIXES: .F90 .F .o .c

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

.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(F90FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(F77FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<



TARGET = h2nc

OBJECTS=\
	m_handle_err.o \
	m_fields_to_plot.o \
	p_h2nc.o 



all: $(TARGET) 

$(TARGET):  $(OBJECTS) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)

install: all
	mkdir -p ../../bin/
	cp $(TARGET) ../../bin


clean:
	rm ./TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(TARGET)



