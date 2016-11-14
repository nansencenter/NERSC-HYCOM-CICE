.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .c

# Include arch-specific makefile macros
# include ../Make.Inc/make.inc
MSCP=../../../../../../..
include $(MSCP)/src/Make.Inc/make.inc

# Library file name
LIB_NERSC   = libhycnersc.a
LIB_CONFMAP = libconfmap.a

CPPFLAGS:=$(CPPFLAGS) 

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
# (called CPPFLAGSF here)
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGSF:=$(subst -,$(subs),$(CPPFLAGS))
else
   CPPFLAGSF:=$(CPPFLAGS)
endif


# Rules for running cpp and updating files in TMP directory
.F90.o:
	$(CF90) -c  $(CPPFLAGSF) $(FFLAGS) $(F90FLG) -o $*.o $<
.F.o:
	$(CF77) -c  $(CPPFLAGSF) $(FFLAGS) $(F77FLG) -o $*.o $<
.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS)  -o $*.o $*.c


# Rules for updating hycom library 
LIB_NERSC_OBJECTS= spherdist.o\
						 wtime.o\
						 sort.o\
						 mod_xc.o\
						 mod_za.o\
						 rotate.o\
						 rotate2.o\
						 machi_c.o\
						 machine.o\
						 mod_parameters.o\
						 mod_year_info.o\
						 mod_hycomfile_io.o \
						 mod_grid.o \
						 mod_spline_calc.o
$(LIB_NERSC) : $(LIB_NERSC)($(LIB_NERSC_OBJECTS))



# Rules for updating confmap library 
LIB_CONFMAP_OBJECTS= mod_confmap.o

$(LIB_CONFMAP) : $(LIB_CONFMAP)($(LIB_CONFMAP_OBJECTS))



lib : $(LIB_HYCOM) $(LIB_NERSC) $(LIB_CONFMAP)

all : lib

install : all
	mkdir -p ../include/
	cp *.mod ../include
	mkdir -p ../lib/
	cp *.a 	 ../lib

clean:
	rm -f *.o *.mod *.a

#dependencies
mod_za.o : mod_xc.o
