#!/usr/bin/env python
import argparse
import modelgrid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile.abfile as abf
import modeltools.cice.io
import numpy
from mpl_toolkits.basemap import Basemap
import logging

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False


def main(grid) :
   griddict=grid.create_datadict_hycom()
   abf.write_regional_grid(griddict)
   abf.write_diag_nc(griddict)
   logger.info("grid shown in grid.png")
   grid.plotgrid(2.).canvas.print_figure("grid.png")

   # Also dump cice grid file
   modeltools.cice.io.write_netcdf_grid(grid,"cice_grid.nc")


def main_proj4(proj4_string,ll_lon,ll_lat,dx,dy,nx,ny):
   # Create grids and write to file
   grid1=modelgrid.Proj4Grid(args.proj4_string,args.ll_lon,args.ll_lat,args.dx,args.dy,args.nx,args.ny)
   main(grid1)

def main_confmap(filename) :
   # Create grids and write to file
   grid1=modelgrid.ConformalGrid.init_from_file(filename)
   main(grid1)




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Prepare HYCOM grid files')
   subparsers = parser.add_subparsers(help='sub-command help')

   parser_confmap = subparsers.add_parser('confmap', help='')
   parser_confmap.add_argument('--filename', type=str, default="grid.info",  help="proj4 string ")
   parser_confmap.set_defaults(subparser_name="confmap")

   parser_proj4   = subparsers.add_parser('proj4', help='')
   parser_proj4.add_argument('proj4_string', help="proj4 string ")
   parser_proj4.add_argument('ll_lon', type=float, help='lower left corner longitude')
   parser_proj4.add_argument('ll_lat', type=float, help='lower left corner latitude')
   parser_proj4.add_argument('dx',     type=int, help='Grid spacing 1st index [m]')
   parser_proj4.add_argument('dy',     type=int, help='Grid spacing 2nd index [m]')
   parser_proj4.add_argument('nx',     type=int, help='Grid dimension 1st index []')
   parser_proj4.add_argument('ny',     type=int, help='Grid dimension 2nd index []')
   parser_proj4.set_defaults(subparser_name="proj4")

   args = parser.parse_args()
   if args.subparser_name == "confmap" :
      main_confmap(args.filename)
   elif args.subparser_name == "proj4" :
      main_proj4(args.proj4_string,args.ll_lon,args.ll_lat,args.dx,args.dy,args.nx,args.ny)
