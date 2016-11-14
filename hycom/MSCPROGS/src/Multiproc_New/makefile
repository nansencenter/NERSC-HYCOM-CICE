VPATH = .:TMP
.SUFFIXES:
.SUFFIXES: .F90 .F .o

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
	cd ./TMP ; $(CF90) -c $(CPPFLAGS)  $(F90FLG) $(FFLAGS) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(F77FLG) $(FFLAGS)  $(INCLUDE) -o $*.o ../$<

TARGET = multiproc

all : $(TARGET)

OBJECTS = mod_pakk.o \
		   m_pakk_field_ids.o\
			p_multiproc.o

install : all 
	mkdir -p ../../bin
	cp $(TARGET) ../../bin

$(TARGET):  $(OBJECTS) 
	cd TMP; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 

clean:
	rm ./TMP/*.o TMP/*.mod $(TARGET)



