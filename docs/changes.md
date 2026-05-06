This page summarizes some of the changes from NERSC HYCOM 2.2.37 to the code in
this version. Most of these changes were introduced to make the code closer to
that of the original code from hycom.

# infile.in

This file is no longer used in this version of the code.

In NERSC HYCOM 2.2.37, infile.in controlled a number of aspects of the
model, such as when to start the model, when to stop, when to output diagnostic,
etc etc. Most of thes have been moved to blkdat.input, or controlled by
preprocessing scripts. Below is a list of the infile.in contents and how they are controlled
now.

|source file     | purpose|
|-------- | -------------|
|version number of run (rungen)  | Not used. Use combination of region name + expt number, or use prefix for archive files etc (can be set in blkdat.input)|
|Reference year, nday1, nday2| Set by expt_preprocess.sh, based on actual dates with iso format | 
|Forcing option climatology option| Moved to preprocessing (script atmo_synoptic.sh) that creates .a and .b files. In-line forcing from netcdf files not supported|
|Surface Temperature relaxation scale |set in thermf.F (def 30 days) |
|Surface Salinity relaxation scale | set in thermf.F (def 30 days) |
|Monthly average | Removed. Set meanfq in blkdat.input for temporal averaging |
|Daily average | Removed. Set meanfq in blkdat.input for temporal averaging |
|lnesto,nesto | Removed. nesting is based on archv files |
|lnesti,nesti | Nesting is based on archv files, output of these is controlled in blkdat.input |
|Tides  | Must be preprocessed by tidal_forcing.sh |
|lgridp  |Removed |
| Diagnostic days | Removed - can be controlled by frequencies in blkdat.input  (diagfq, dsurfq)

# Daily averages

The daily average module is no longer used. In stead, you should use the mean
diagnostics option in blkdat.input (Look for meanfq). Further postprocessing (ex
ensemble staticsics) must be done using separate programs.

# Weekly averages

The weekly average module is no longer used. In stead, you should use the mean
diagnostics option in blkdat.input (Look for meanfq). Further postprocessing (ex
ensemble staticsics) must be done using separate programs.


# Nesting

Nesting is now done using standard hycom routines. The outer nesting region must
be run first to create archive files at fixed intervals, these must then be
interpolated horizontally (and perhaps also vertically) to create nesting
conditions to the inner models.  More info here (TODO)[TODO]

# Tides 

Tides are now implemented by adding barotropic pressure terms and velocities to
existing nesting fields. 
More info here (TODO)[TODO]

# mod_year_info

This module used to handle calendars in the NERSC hycom version, but has been removed
in this version. 

# Naming of outout files

Previous version of NERSC hycom used hardcoded file names with a prefix set in
infile.in and possibly an ensemble member provided as input argument to the
hycom executable. This version uses neither of these. In stead, it is possible
to set a file prefix directly in blkdat.input.
