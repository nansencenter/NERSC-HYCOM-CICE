# makefile includes for hexagon, gnu compilers
# Standard compilers and linkers
CF90= ftn
CF77= ftn
CC= cc
LD= ftn

# CPP, used internally by compilers - note that some compilers
# need special options for this, consult compiler manual if in trouble
CPP= /usr/bin/cpp -P -traditional

#NCAR graphics compiler wrappers
NCARGCF90=ncargf90
NCARGCF77=ncargf77
NCARGCC=ncargcc
NCARGLD=ncargf90

# Flags for fixed/free format
F90FLG= -ffree-form
F77FLG= -ffixed-form

# Compiler flags, real8 version
FFLAGS= -O2 -fconvert=big-endian -g -fdefault-real-8
CFLAGS= -O2

# Compiler flags, real4 version (needed for NCARG routines)
FFLAGSR4= -O2 -fconvert=big-endian -g -fdefault-real-4
CFLAGSR4= -O2

# Linker flags
LINKFLAGS= $(FFLAGS)  

#Netcdf, FFTW and lapack Libraries. Not needed on hexagon
INCLUDE= 
LIBS = 

# Some fortran compilers have iargc as built in, 
# others as library routine
CPPFLAGS=-UIARGC -DFFTW -DLAPACK

