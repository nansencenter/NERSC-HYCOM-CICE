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
	cd ./TMP ; $(CF90) $(CPPFLAGS) $(F90FLG) $(FFLAGS) $(INCLUDE) -c -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) $(CPPFLAGS) $(F77FLG) $(FFLAGS) $(INCLUDE) -c -o $*.o ../$<



TARGET = tec_nc_conv
scripts= m2nc m2t

#These must be in correct order
OBJECTS=\
	mod_types.o\
	m_handle_err.o \
	m_tecplot_header.o\
	m_fields_to_plot.o  \
	m_tecplot_zoneinfo.o\
	m_tecplot_dump.o\
	m_tecplot_dump_rot.o \
	tecconv.o\
	lon_lat_extrap.o



all: $(TARGET) 

install: all
	mkdir -p ../../bin/
	cp $(scripts) ../../bin
	cp $(TARGET) ../../bin

$(TARGET):  $(OBJECTS) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)

clean:
	rm TMP/*.o TMP/*.mod $(TARGET)



