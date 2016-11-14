
##############################################################################
#       ISLSCP2 Land/Water and Land/Sea Masks and other Ancillary Data       #
#                at 0.25, 0.5 and 1 degree Spatial Resolutions               #
##############################################################################

The files in this directory were produced by the ISLSCP2 staff based on the 30 
arcsecond (~1km) Earth Observing System (EOS) Digital Elevation Model (DEM) 
provided by Dr. Thomas Logan from NASA/Jet Propulsion Laboratory (JPL). These 
are global, coarse resolution (1/4, 1/2 and 1 degree spatial resolutions) binary 
land/water and land/sea masks (i.e. land or water/sea), and binary masks of 
inland water bodies; and files with the fraction of land, water, sea or inland water 
within each cell.

Many global studies and/or models require a priori knowledge of the spatial 
distribution of land masses and water bodies on the surface of the Earth. In 
support of the International Satellite Land Surface Climatology Project (ISLSCP) 
Initiative II user community, we have created a series of land/water and 
land/sea masks along with other potentially useful ancillary data sets. The 
original intent was to provide these data as separate data sets so that the user 
would have some flexibility in applying any mask to the various data sets of the 
ISLSCP2 data collection. Based on recommendations from the ISLSCP2 Science 
Working Group and the modeling community, the ISLSCP2 land/water and/or land/sea 
masks have been applied to all ISLSCP2 data sets that require clear land/water 
boundaries. Users should use these files as a resource should any questions 
arise with land/sea or land/water boundaries. These data were used to make the 
entire ISLSCP2 data collection internally consistent based on the needs of the 
modeling community. However, users can use the percentages of land/water/sea to 
create different masks based on percentages that may better address their 
particular needs.

The ISLSCP2 masks are based on the best available 30 arcsecond (~1km) land/water 
masks produced at JPL in support of the EOS AM-1 satellite platform. The JPL 
masks are based in turn on vector data from the World Vector Shoreline (WVS) 
database (Solluri and Woodson 1990) for coastlines, and the Digital Chart of the 
World (DCW) (Danko 1992) for inland water bodies. The WVS is based on older 
1:250,000 scale Defense Mapping Agency JOG (Joint Operational Graphics) charts 
produced from a variety of sources. DCW is based on a photogrammetric consistent 
mapping at a scale of 1:1,000,000. The ISLSCP2 staff have aggregated the 1km 
data to spatial resolutions of 1/4, 1/2 and 1 degree, in the process generating 
layers with the percentage of land/water or land/sea cells within the coarse 
resolution cells. An arbitrary threshold of greater than or equal to 50% water 
was used to designate water dominated cells and to produce binary land/water or 
land/sea masks for the data collection.

Please acknowledge the EOS DEM Science Working Group and the EOS ASTER 
Instrument Project at JPL for provision of the 30 arcsecond EOS AM1 DEM and the 
ISLSCP2 data collection when these data are used. 

Any problems encountered with the data set should be reported to Eric Brown de 
Colstoun at ericbdc@ltpmail.gsfc.nasa.gov

References:
----------
Danko, D.M. (1992). The Digital Chart of the World Project, Photogrammetric 
Engineering and Remote Sensing, 58:1125-1128.

Soluri, E.A. and V.A. Woodson (1990) World Vector Shoreline. International 
Hydrographic Review LXVII(1).

##############################################################################

Directory Listing
=================
The data sets in this directory are provided at three spatial resolutions of 
0.25, 0.5 and 1 degrees in latitude and longitude. The following files are 
included:

Standard Land/Water Masks
Contained in the file "land_water_masks_xdeg.zip"
-------------------------------------------------
  land_water_mask_xx.asc: Global gridded binary land/water masks, with values of 
  0 and 1 (Water=0; Land=1). xx can be qd, hd, and 1d, denoting a spatial 
  resolution of 1/4, 1/2 and 1 degree, respectively.

  land_percent_xx.asc: Percentage of original 30 arcsecond land cells located 
  within coarser resolution cells, from 0 to 100. xx is the resolution as above.
  These files were used to construct the binary land/water masks.

  water_percent_xx.asc: Percentage of original 30 arcsecond water cells (including 
  inland water) located within coarser resolution cells, from 0 to 100. xx is the 
  resolution as above. These files were used to construct the binary land/water masks.


Special Land/Water Masks (land=-88, water=-99)
Contained in the file "land_water_masks-99_xdeg.zip"
----------------------------------------------------
  land_water_mask-99_xx.asc: Global gridded binary land/water masks, with values
  of -88 and -99 (Water=-99; Land=-88). xx is the resolution as above. These
  files were created from the "land_water_mask_xx.asc" files listed above. These files
  were used as the bases for several derived map products (point data was plotted
  on top of these maps), so they are included here as a baseline comparison product.


Special Land/Ocean Masks (no inland water bodies)
Contained in the file "land_ocean_masks_xdeg.zip"
-------------------------------------------------
  land_ocean_mask_xx.asc: Global gridded binary land/ocean masks, with values of 0 
  and 1 (Ocean=0; Land=1). xx can be qd, hd, and 1d denoting a spatial resolution of 
  1/4, 1/2 and 1 degree, respectively. NOTE that these are masks of only land 
  and oceans and do not contain inland water bodies.

  land_percent2_xx.asc: Percentage of original 30 arcsecond land cells located 
  within coarser resolution cells, from 0 to 100. xx is the resolution as above.
  These files were used to construct the binary land/ocean masks. NOTE that these
  files are different from the "land_percent_xx.asc" files because inland water
  bodies are not considered.

  ocean_percent2_xx.asc: Percentage of original 30 arcsecond sea or ocean cells 
  located within coarser resolution cells, from 0 to 100. xx is the resolution as 
  above. These files were used to construct the binary land/ocean masks. NOTE that
  these files do not contain inland water bodies and are different from the
  "water_percent_xx.asc" files.


Special Inland Water Masks (only inland water bodies)
Contained in the file "inland_water_masks_xdeg.zip"
---------------------------------------------------
  inland_water_mask_xx.asc: Global gridded binary inland water bodies masks, with 
  values of 0 and 1 (Inland Water=1; Other=0). xx can be qd, hd, and 1d denoting 
  a spatial resolution of 1/4, 1/2 and 1 degree, respectively. NOTE that these are
  masks of only inland water bodies, no land or oceans are included.

  other_percent_xx.asc: Percentage of original 30 arcsecond non-inland water cells 
  (i.e. land and oceans) located within coarser resolution cells, from 0 to 100. 
  xx is the resolution as above. These files were used to construct the binary
  inland water bodies masks.

  inland_water_percent_xx.asc: Percentage of original 30 arcsecond inland water 
  cells located within coarser resolution cells, from 0 to 100. xx is the resolution
  as above. These files were used to construct the binary inland water bodies masks.


##############################################################################

ASCII File Format
------------------
All of the files in the ISLSCP Initiative II data collection are in the ASCII, 
or text format. The file format consists of numerical fields of varying length, 
which are delimited by a single space and arranged in columns and rows. The 
values for the binary land/water or land/sea masks are written as the integers
0 and 1. All values in the land/water/sea fraction files are written as
integers from 0 to 100.

The files at different spatial resolutions each contain the following numbers of 
columns and rows:

   One degree (1d): 360 columns by 180 rows 
   1/2 degree (hd): 720 columns by 360 rows 
   1/4 degree (qd): 1440 columns by 720 rows 

All files are gridded to a common equal-angle lat/long grid, where the 
coordinates of the upper left corner of the files are located at 180 degrees W, 
90 degrees N and the lower right corner coordinates are located at 180 degrees E,
90 degrees S. Data in the files are ordered from North to South and from West to
East beginning at 180 degrees West and 90 degrees North.

The data files are PKZip compressed. On UNIX, they can be decompressed using the 
"unzip -a" command. On Windows, they can be decompressed using WinZip or other 
PKZip software. On the Macintosh, they can be decompressed using Stuffit 
Expander.

################################################################################
