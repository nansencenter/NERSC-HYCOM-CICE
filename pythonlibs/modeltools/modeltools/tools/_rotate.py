import numpy



# This class is the same as function below, but splits into preprocessing and actual rotation. Faster...
class rotateVector(object):

   def __init__(self,mlon,mlat) :
      # Rotates a vector field from a grid in geographical coordinates
      # into or from the grid defined by the lat,lon in the input variables.
      # C-grid is assumed.                     
      #  
      # Input: mlat, mlon: position in scalar point
      # ----------------------------------------------------------------------
      radian = 180./numpy.pi
      radinv=1./radian
      Itest=50
      Jtest=50
      # ----------------------------------------------------------
      # Assumes that all parameters are provided in P point,
      # perform the rotation rquired in curvlinear grid.
      # -------------------------------------------------------   
      # Rotation angle in P-point 
      dlon=numpy.zeros(mlon.shape)
      dlat=numpy.zeros(mlon.shape)
      dlon[:,1:-1]=mlon[:,2:] - mlon[:,:-2]  # Forward azimuth
      dlat[:,1:-1]=mlat[:,2:] - mlat[:,:-2]  # Forward azimuth
      dlon[:,-1 ]=dlon[:,-2]
      dlat[:,-1 ]=dlat[:,-2]
      dlon[:,0 ]=dlon[:,1]
      dlat[:,0 ]=dlat[:,1]
      dlon = numpy.mod(dlon+360+180.,360.)-180 # Make sure in range (-180, 180)
      weight_lat=numpy.cos(radinv*mlat)
      #
      # Angle i-dir
      self._theta_up = numpy.arctan2(dlat,dlon*weight_lat)   
      #
      # Rotation angle in P-point 
      dlon=numpy.zeros(mlon.shape)
      dlat=numpy.zeros(mlon.shape)
      dlon[1:-1,:]=mlon[2:,:] - mlon[:-2,:]  # Forward azimuth
      dlat[1:-1,:]=mlat[2:,:] - mlat[:-2,:]  # Forward azimuth
      dlon[-1,:]=dlon[-2,:]
      dlat[-1,:]=dlat[-2,:]
      dlon[0,:]=dlon[1,:]
      dlat[0,:]=dlat[1,:]
      dlon = numpy.mod(dlon+360.+180.,360.)-180 # Make sure in range (-180, 180)
      weight_lat=numpy.cos(radinv*mlat)
      #
      # Angle j-dir
      self._theta_vp = numpy.arctan2(dlat,dlon*weight_lat)

   def rotate(self,ud,vd) :
      # Rotated velocities
      urot= ud*numpy.cos(self._theta_up)+ vd*numpy.sin(self._theta_up)
      vrot= ud*numpy.cos(self._theta_vp)+ vd*numpy.sin(self._theta_vp)
      return urot,vrot




def rotate_vector(ud,vd,mlon,mlat) :
   # Rotates a vector field from a grid in geographical coordinates
   # into or from the grid defined by the lat,lon in the input variables.
   # C-grid is assumed.                     
   #  
   # Input: mlat, mlon: position in scalar point
   #        nx,ny:     dimension of the model grid.
   #        ud,vd:  Unrotated vector components,  where ud is the EW
   #                component and vd is the NS component 
   #        dir:    l2m (latlon to general)
   #                m2l (general to latlon)
   #
   # Output: ud,vd: Rotated vector components, where  ud is along the
   #                i-axis and vd is along the j-axis.
   # ----------------------------------------------------------------------
   radian = 180./numpy.pi
   radinv=1./radian
   Itest=50
   Jtest=50

   # ----------------------------------------------------------
   # Assumes that all parameters are provided in P point,
   # perform the rotation rquired in curvlinear grid.
   # -------------------------------------------------------   

   # Rotation angle in P-point 
   dlon=numpy.zeros(mlon.shape)
   dlat=numpy.zeros(mlon.shape)
   dlon[:,1:-1]=mlon[:,2:] - mlon[:,:-2]  # Forward azimuth
   dlat[:,1:-1]=mlat[:,2:] - mlat[:,:-2]  # Forward azimuth
   dlon[:,-1 ]=dlon[:,-2]
   dlat[:,-1 ]=dlat[:,-2]
   dlon[:,0 ]=dlon[:,1]
   dlat[:,0 ]=dlat[:,1]
   dlon = numpy.mod(dlon+360+180.,360.)-180 # Make sure in range (-180, 180)
   weight_lat=numpy.cos(radinv*mlat)

   
   # Angle i-dir
   theta_up = numpy.arctan2(dlat,dlon*weight_lat)   

   # Rotation angle in P-point 
   dlon=numpy.zeros(mlon.shape)
   dlat=numpy.zeros(mlon.shape)
   dlon[1:-1,:]=mlon[2:,:] - mlon[:-2,:]  # Forward azimuth
   dlat[1:-1,:]=mlat[2:,:] - mlat[:-2,:]  # Forward azimuth
   dlon[-1,:]=dlon[-2,:]
   dlat[-1,:]=dlat[-2,:]
   dlon[0,:]=dlon[1,:]
   dlat[0,:]=dlat[1,:]
   dlon = numpy.mod(dlon+360.+180.,360.)-180 # Make sure in range (-180, 180)
   weight_lat=numpy.cos(radinv*mlat)

   # Angle j-dir
   theta_vp = numpy.arctan2(dlat,dlon*weight_lat)

   # Rotated velocities
   urot= ud*numpy.cos(theta_up)+ vd*numpy.sin(theta_up)
   vrot= ud*numpy.cos(theta_vp)+ vd*numpy.sin(theta_vp)

   #print mlat.min(),mlat.max()
   #print "mlat.shape :",mlat.shape
   #print "weight_lat :",weight_lat[Itest,Jtest]
   #print "dlon       :",dlon[Itest,Jtest]
   #print "dlat       :",dlat[Itest,Jtest]
   #print "mlat       :",mlat[Itest,Jtest]
   #print "theta1R    :",theta_up[Itest,Jtest]*radian
   #print "theta2     :",theta_vp[Itest,Jtest]*radian

   return urot,vrot

