# ######################################################################
#  Makefile FES prediction software.
# 
#  File      : Makefile
#  Developer : CLS
#  Version   : 1.5
#  Date      : 22 March 2007
#  
#  This software have been tested on Linux platform
#  It is provided without any warranty.
# 
#  For bug reports, please contact :
#  ---------------------------------
#  Fabien LEFEVRE 
# 
#  CLS
#  http://www.cls.fr
#  Direction Océanographie Spatiale
#  8-10, rue Hermès - Parc Technologique du Canal
#  31520 Ramonville Saint-Agne cedex - France
#  Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
#  e-mail: Fabien.Lefevre@cls.fr
# 
#  NOTE: This software is based on the former versions
#        developed by CNRS (Jean-Marc MOLINES and
#        Florent LYARD among others)
#
#
#   IK: version running on Fimm 
#  KAL: Mild cleanup of Inti's version - moved fes routines and library 
#       Out of this directory. This is cleaner when the database is updated
# #####################################################################

# --------------------------------------------------------------------
# TO BE MODIFIED IF NEEDED !!!
# 
#CFLAGS=-march=pentium4 -mfpmath=sse -msse2 -mmmx -O3 -pipe -fomit-frame-pointer -Wall -Wpointer-arith -Wshadow -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Werror
#FFLAGS=-O
# CFLAGS=-g
 CFLAGS=-O
# FFLAGS=-g
#
#-- If operating system is Linux
# CC		= gcc
 FC		= pgf77
#
#-- If operating system is Sun Solaris
#CC 		= cc
#FC		= f77
#
# End of TO BE MODIFIED IF NEEDED !!!
# --------------------------------------------------------------------

#INCLUDE = -I/home/fimm/nersc/intissar/Progs/Tide_FES/include
#LIB_DIR =/home/fimm/nersc/intissar/Progs/Tide_FES/lib

# Knut - use relative paths - Now to FES2004 in separate dir
# NB: must use gnu-linked netcdf ( or change compilers, use Make.Inc?)
NETCDF_LIB=-L/local/netcdf-3.6.2-gnu/lib -lnetcdf 
INCLUDE = -I../../../include  -I/local/netcdf-3.6.2-gnu/include 
LIB_DIR = ../../../lib
LIB_FES = $(LIB_DIR)/libfes.a
LDFLAGS = -L$(LIB_DIR)


.SUFFIXES: .o .f .c

.c.a:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $*.o
	$(AR) r $@ $*.o
	$(RM) $*.o

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $<

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE) $<


BINARIES = icebfes fes2mod fes8const fes2nc fes2hycnc fes2nc_transposed

#
# all
#
all : $(BINARIES)

#
# clean
#
clean :
	$(RM) *.o $(BINARIES) fes2mod fes2nc fes8const icebfes

#
#  fes2mod
#
#
fes2mod : fes2mod.o fesCore8Const.o
	$(CC) -o $@ $(CFLAGS) fes2mod.o fesCore8Const.o $(LDFLAGS) -lfes -lm  

#  fes2nc
#
fes2nc : fes2nc.o  fesCore8Const.o
	$(CC) -o $@ $(CFLAGS) fes2nc.o fesCore8Const.o $(LDFLAGS) -lfes -lm  $(NETCDF_LIB)
#
#  fes2nc_transposed
#
fes2nc_transposed : fes2nc_transposed.o  fesCore8Const.o
	$(CC) -o $@ $(CFLAGS) fes2nc_transposed.o fesCore8Const.o $(LDFLAGS) -lfes -lm  $(NETCDF_LIB)
#
#  fes2ncts
#
fes2hycnc : fes2hycnc.o  fesCore8Const.o
	$(CC) -o $@ $(CFLAGS) fes2hycnc.o $(LDFLAGS) -lfes -lm  $(NETCDF_LIB)


#  fes8const
#
fes8const : fes8const.o  fesCore8Const.o
	$(CC) -o $@ $(CFLAGS) fes8const.o  fesCore8Const.o  $(LDFLAGS) -lfes -lm	

#
# icebFes 
#
icebfes : icebfes.o 
	$(CC) -o $@ $(CFLAGS) icebfes.o $(LDFLAGS) -lfes -lm	



install : all
	mkdir -p ../../../bin_setup/
	cp $(BINARIES) ../../../bin_setup/
