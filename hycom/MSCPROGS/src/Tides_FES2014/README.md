

The main routine here for nersc-hycom-cice is the fes2hycom.c routine.
This routine will read FES2014 data, calculate tides and write tidal elevation 
to a netcdf file. so, output is tidal timeseries, not constituents....

Run like this
   fes2hycom h 1.0 1.25 1.5 1.75 ....

Input is time in days since 1950-01-01 00:00:00 UTC. First letter indicates if you wish to calculate elevations(h), u-component (u) or  v-component (v)
of current. Output is a netcdf file.

fes2hycom reads the time info, and it will also read  regional.grid and regional.depth 
files to get the locations where tidal input is needed (these need to be present in your
working directory). 

In addition you will need to set the path to the FES 2014  data. This path can be set 
as an environment variable, something like this:

    export FES2014_PATH=/path/to/FES2014/data

fes2hycom.c will in general give you some information of what went wrong if there are errors. 


# FES library

**NB: Dont use the FES library that is accompanying MSCPROGS (in src_others). It uses an older version of the FES C API.**

You will need to have the FES 2014 C API installed to use these routines. This can be found on github at
https://bitbucket.org/cnes_aviso/fes

Download and install from  there, or have an IT guy install it. I recommend to use the gnu compiler to compile, as this 
plays best with autoconf etc...


# FES memory use

fes2hycom.c is set up to  use FES_MEM as reading method when fesNew is called. This loads
all tidal constituenbts into memory, which may be problematic if you have little memory on your machine. If 
you run out of memory, try to call fes_new with FES_IO in stead of FES_MEM. This will take a long time for 
a large grid, but uses less memory. 

In practice, loading all tidal constituents when using FES_MEM uses around 4 GB, which should not be a big problem on a modern machine.

Note that FES_IO takes a long time for the first time you call it. Subsequent time steps are much faster due to a caching routine in FES C api.



# Compilation

Compile with gnu c compiler (gcc, or cc with PrgEnv-gnu on hexagon). See makefile for more info. 



