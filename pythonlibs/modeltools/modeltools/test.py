import unittest
import modeltools.grid

class GridTest(unittest.TestCase):
    p4s="+proj=stere +lon_0=0 +lat_0=60 +lat_ts=60 +ellps=sphere"
    ll_lon=-20
    ll_lat=60
    dx=4000
    dy=4000
    nx=500
    ny=300

    def test_grid_shape(self):
       grid1=modeltools.grid.Proj4Grid(self.p4s,self.ll_lon,self.ll_lat,self.dx,self.dy,self.nx,self.ny)
       plon,plat=grid1.pgrid()
       self.failUnlessEqual(plon.shape[1],self.nx)
       self.failUnlessEqual(plon.shape[0],self.ny)
       self.failUnlessEqual(plat.shape[1],self.nx)
       self.failUnlessEqual(plat.shape[0],self.ny)


    def test_ll_inversion(self) :
       grid1=modeltools.grid.Proj4Grid(self.p4s,self.ll_lon,self.ll_lat,self.dx,self.dy,self.nx,self.ny)
       plon,plat=grid1.pgrid()
       self.assertAlmostEqual(plon[0,0],self.ll_lon,places=7)
       self.assertAlmostEqual(plat[0,0],self.ll_lat,places=7)

       
if __name__ == "__main__" :
   unittest.main()
       
