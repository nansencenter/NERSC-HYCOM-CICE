

The main routine here for nersc-hycom-cice is the fes2004hycom.c routine.
This routine will read FES2004 data, calculate tides and write tidal elevation 
to a netcdf file. so, output is tidal timeseries, not constituents....

Run like this
   fes2004hycom 1.0 1.25 1.5 1.75 ....

Input is time in days since 1950-01-01 00:00:00 UTC.

fes2004hycom reads the time info, and it will also read  regional.grid and regional.depth 
files to get the locations where tidal input is needed (these need to be present in your
working directory). 

In addition you will need to set the path to the FES 2004  data. This path can be set 
as an environment variable, something like this:

    export FES2004_PATH=/path/to/FES2004/data

fes2004hycom.c will in general give you some information of what went wrong if there are errors. 


# Other routines

The remaining routines in this directopry were mainly used to calculate tidal consituents
for use in older hycom version. These are not needed when running the new hycom-cice version, as
the new version uses tidal elevations directly. 

Note that if you try to run these routines you can run into problems because of stacksize limityatiuons. If
that happens, try to run this command before running 

   ulimit -s unlimited

# FES memory use

fes2004hycom.c is set up to  use FES_MEM as reading method when fesNew is called. This loads
all tidal constituenbts into memory, which may be problematic if you have little memory on your machine. If 
you run out of memory, try to call fes_new with FES_IO in stead of FES_MEM. This will take a long time for 
a large grid, but uses less memory.

In practice, loading all tidal constituents uses around 500MB - 1 GB, which should not be a big problem on a modern machine.


# Compilation

Compile with pgi c compiler (cc with PrgEnv-pgi on hexagon).  gnu compilers should also work, but wil give lots of errors. It is possible to suppress these with compiler flags.




# Old  documentation. 

Below is some old documentation for FES2004. It may not be 100% updated.




This routine contains routines for producing tidal forcing of
the hycom model, and for visualizing fes data (mainly fes2nc).

Make sure the environment variable "FES_PATH" is set when you 
run the routines. It points to the directory containing the FES 
data. (Relative to this dir it is ../FES2004/data).

The main routines are:

fes2mod : Produces input files for forcing hycom
fes2nc  : Puts FES data into a netcdf file for visualization

See ../Inti_README.txt for more info on running fes2mod

Also, see Inti_notes.doc in this directory for more info, altough
some of that isnt 100% correct in this setup (some refer to routines
in ../FES2004/src)


for fes2nc:
Make sure that netcdf libraries and compilers are consistent. 
Use gnu compilers and libraries if present

If it doen not compile, it may be that libfes.a is not yet compiled, try this:
cd hycom/MSCPROGS/src_others/FES2004/src/
make clean
make install
cd ~/hycom/MSCPROGS/src/Tides_FES/nersc_src
make clean
make install
