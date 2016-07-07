#
# Makefile for archv2nc
#
# --- usage: make archv2nc ARCH=sun >& Make_archv2nc
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#
# --- assumes netcdf module and library are in this directory
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1
.c:
	@echo "Must have an explicit rule for" $*
	@exit 1

include ../../config/$(ARCH)_setup

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_ncdf.com

clean:
	/bin/rm -f *.o *.a *.mod *2nc

MODS =  date_sub.o globals.o 
OBJS =  ncrange.o wrtNcdf.o archv2nc.o raw.o locate.o
SRCS =	archv2nc.f wrtNcdf.f ncrange.f locate.f raw.f globals.F date_sub.F  

#
# --- executables.
#

archv2nc:	$(MODS) $(OBJS)
	$(LD)  $(LDFLAGS) -o archv2nc  $(MODS) $(OBJS) $(EXTRALIBS) libnetcdf.a $(EXTRANCDF)

#
# --- explicit dependencies for each .[fF] file.
#

archv2nc.o: archv2nc.f globals.o date_sub.o wrtNcdf.o raw.o locate.o
wrtNcdf.o: wrtNcdf.f ncrange.o globals.o
ncrange.o: ncrange.f
globals.o: globals.F
data_sub.o: date_sub.F
raw.o: raw.f
