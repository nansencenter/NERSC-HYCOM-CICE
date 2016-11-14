# Procedure for compiling programs in hycom_2.2.72_ALL.

## Prerequisites

You will need to have cshell installed to runn the make scripts

## Instructions

* Go into directory hycom_2.2.72_ALL

* Look inside the config directory to find a configuration for your system. gfortran_setup will work on machines with a gfortran compiler, amd64 will work on machines with a portland compiler

* Most routines dont rely on external libs, but some netcdf conversion routines will need to know where the netcdf libraries are. If you need these, make sure they are set in the config file

* When you have found a setup, edit Make_all.src. There you should set ARCH to the name of the config file minus the "_setup" part. For instance, to use gfortran_setup, it should read

    setenv ARCH gfortran
    
* Also note that there is a make script inside the hycom_2.2.72_ALL/bin/directory that has to be edited by hand. This one is a bit messy, but look for a suitable setup for your architecture in one of the case statements...
    
* After you edited Make_all.src, run the make scripts like this:

    csh ./Make_all.com
    
* When the compile script runs, you will get feedback about the progress. Programs will be created in subdirs archive force meanstd plot relax sample subregion topo.

* There will also be several useful routines in the bin directory created by the compile script there,