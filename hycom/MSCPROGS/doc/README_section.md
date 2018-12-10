# Extract section data and tranports

The **Section/** subdirectory contain routines for extracting data from sections to NetCDF files, and routines for calculating transport estimates from hycom files. There are several programs in this directory, but their functionality is most easily accessed using three different scripts; they are **m2section**, **m2transport** and **m2transport2**.

All the scripts in this directory will first call a routine **section_intersect**, whose purpose is to set up the grid points which are needed to calculate transports and to interpolate section details. The sections are specifed in the file [sections.in](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/master/hycom/MSCPROGS/Input/sections.in). Once a section is read, the grid points along a section are calculated.
