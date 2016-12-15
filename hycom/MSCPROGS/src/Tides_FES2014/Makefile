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
#  KAL 201612. Modified for FES 2014. Needs external library, see  https://knutal@bitbucket.org/cnes_aviso/fes.git
#              FES_DIR below must point to the root of the external library
#              For now this program only works with the gcc compiler. Mainly because the 
#              version of the fes library I used was compiled with gcc..
# #####################################################################

# --------------------------------------------------------------------
# TO BE MODIFIED IF NEEDED !!!
# 
# on hexagon, cc works when PrgEnv-gnu is loaded. Otherwise use gcc (but you may have to modify stoff below)
# CC = gcc
CC		= cc

#CFLAGS=-march=pentium4 -mfpmath=sse -msse2 -mmmx -O3 -pipe -fomit-frame-pointer -Wall -Wpointer-arith -Wshadow -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Werror
# CFLAGS=-g
CFLAGS= -g    -Wall #-Wpointer-arith -Wshadow -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Werror



# Set these appropriately for your systemn (if needed)
#NETCDF_LIB=-L/local/netcdf-3.6.2-gnu/lib -lnetcdf 
#NETCDF_INC=-I/local/netcdf-3.6.2-gnu/include
NETCDF_LIB=-lnetcdf 
NETCDF_INC=

# Location of fes 
FES_DIR=/home/nersc/knutali/opt/fes-2.7.4/

# End of TO BE MODIFIED IF NEEDED !!!
# --------------------------------------------------------------------

# Normally  need to modify these 
INCLUDE = -I${FES_DIR}/include/ ${NETCDF_INC}
LIBS = -L${FES_DIR}/lib/ -lfes ${NETCDF_LIB}
LDFLAGS = $(LIBS) 


.SUFFIXES: .o .f .c

.c.a:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $*.o
	$(AR) r $@ $*.o
	$(RM) $*.o

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $<

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE) $<


BINARIES = icebfes fes2014hycom 

#
# all
#
all : $(BINARIES)

#
# clean
#
clean :
	$(RM) *.o $(BINARIES) 
#
# fes2hycnc
#
fes2014hycom : fes2014hycom.o 
	$(CC) -o $@ $(CFLAGS) fes2014hycom.o $(LDFLAGS)


#
# icebFes 
#
icebfes : icebfes.o 
	$(CC) -o $@ $(CFLAGS) icebfes.o $(LDFLAGS)


install : all
	mkdir -p ../../bin_setup/
	cp fes2014hycom ../../bin_setup/
