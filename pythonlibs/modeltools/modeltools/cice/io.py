import numpy
import scipy.io.netcdf

def write_netcdf_grid(grid,filename) :
   nc=scipy.io.netcdf.netcdf_file(filename,"w")
   nc.createDimension("x",grid.Nx)
   nc.createDimension("y",grid.Ny)

   nc.createVariable("ulat","double",("y","x",))
   nc.createVariable("ulon","double",("y","x",))
   nc.createVariable("hte","double",("y","x",))
   nc.createVariable("htn","double",("y","x",))
   nc.createVariable("hus","double",("y","x",))
   nc.createVariable("huw","double",("y","x",))
   nc.createVariable("angle","double",("y","x",))

   nc.variables["ulat"].units="radian"
   nc.variables["ulon"].units="radian"
   nc.variables["hte"].units="cm"
   nc.variables["htn"].units="cm"
   nc.variables["hus"].units="cm"
   nc.variables["huw"].units="cm"
   nc.variables["angle"].units="radian"

   # Hycom q point is lower left corner of t-cell.
   # Cice  u point is in pupper left corner of t-cell
   # So we have to fetch the extended grid here

   ulon,ulat=grid.cice_ugrid()
   #print ulon.shape
   #print grid.Nx
   hte=grid.cice_hte()
   htn=grid.cice_htn()
   hus=grid.cice_hus()
   huw=grid.cice_huw()

   # from CICE doc:  ANGLE = angle between local x direction and true east
   azi=grid.p_azimuth() # Forward azimuth of x direction (bearing)
   ang=90. - azi         # Rel. east


   nc.variables["ulat"][:]=numpy.radians(ulat)
   nc.variables["ulon"][:]=numpy.radians(ulon)
   nc.variables["hte"][:]=hte*100. # m -> cm
   nc.variables["htn"][:]=htn*100. # m -> cm
   nc.variables["hus"][:]=hus*100. # m -> cm
   nc.variables["huw"][:]=huw*100. # m -> cm
   nc.variables["angle"][:]=numpy.radians(ang)

   nc.close()



def write_netcdf_kmt(field,filename) :
   nc=scipy.io.netcdf.netcdf_file(filename,"w")
   nc.createDimension("x",field.shape[1])
   nc.createDimension("y",field.shape[0])
   nc.createVariable("kmt","float",("y","x",))
   nc.variables["kmt"].units="1"
   nc.variables["kmt"][:]=field[:]
   nc.close()

