VPATH= .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


include ../Make.Inc/make.inc


# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc  $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/  $(INCLUDE)


.f.o:
	cd TMP; $(CF77) -c  $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$*.f

OBJS=	\
	fillup.o \
	rNODC.o \
	tecle2.o \
	intpol.o\
	reddat.o\
	levini.o \


TARGET=old_levini

$(TARGET):	$(OBJS)
	cd TMP; $(CF77)  -o ../$(TARGET) $(LINKFLAGS) $(OBJS) $(LIBS)

all : $(TARGET)

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET)   ../../bin_setup/

clean:
	rm TMP/*.o $(TARGET)
