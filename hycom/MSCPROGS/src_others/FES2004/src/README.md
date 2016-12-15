Note that the library and include files are no longer copied to MSCPROGS/lib and MSCPROGS/include. This is because FES2014 support is now included, and the FES2014 C
API has changed relative to the FES2004 API. 

The FES2004 include and lib are now copied to ../ . The nersc FES 2004 routines reflect this

# Compilation.

gnu compiler and pgi works. pgi recommended because nersc routines that use FES 2004 give lots of errors with the gnu compiler ...

make lib ; make install should do the trick, unless you need to create the bin files. see makefile for more info, or Readme.txt in parent directory.

