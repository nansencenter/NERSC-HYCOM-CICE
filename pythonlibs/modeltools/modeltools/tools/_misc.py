import logging
import numpy
import pyproj

_default_threshold=-5.
_default_S=0.25

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


def shapiro_filter(w,threshold=_default_threshold,S=_default_S) :
   logger.info("Applying shapiro filter (1 pass, S=%.3f)"%S)
   myw = numpy.copy(w)
   # TODO: Account for periodic grids

   # x - direction
   wm1=numpy.copy(myw[1:-1,1:-3])
   w0 =numpy.copy(myw[1:-1,2:-2])
   wp1=numpy.copy(myw[1:-1,3:-1])
   fwm1=numpy.where(wm1<threshold,1,0)
   fw  =numpy.where(w0 <threshold,1,0)
   fwp1=numpy.where(wp1<threshold,1,0)
   I = numpy.where(fwm1+fw+fwp1==3)
   if I[0].size > 0 :
      w0[I] = w0[I] + 0.5*S * (wm1[I]+wp1[I]-2*w0[I])
   myw[1:-1,2:-2] = w0


   # x - direction
   wm1=numpy.copy(myw[1:-3,1:-1])
   w0 =numpy.copy(myw[2:-2,1:-1])
   wp1=numpy.copy(myw[3:-1,1:-1])
   fwm1=numpy.where(wm1<threshold,1,0)
   fw  =numpy.where(w0 <threshold,1,0)
   fwp1=numpy.where(wp1<threshold,1,0)
   I = numpy.where(fwm1+fw+fwp1==3)
   if I[0].size > 0 :
      w0[I] = w0[I] + 0.5*S * (wm1[I]+wp1[I]-2*w0[I])
   myw[2:-2,1:-1] = w0

   return myw




def remove_one_neighbour_cells(inv,threshold=_default_threshold,num_neighbour=1) :
   """ Sets cells with values above "threshold" to threshold value if it has less  
       than or equal to num_neighbour neighbours """

   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[0,:] = threshold
   v[:,-1] = threshold
   v[:,0] = threshold
   
   # Test
   #v[100,100]=threshold+1
   #v[100, 99]=threshold-1
   #v[100,101]=threshold-1
   #v[ 99,100]=threshold-1
   #v[101,100]=threshold-1

   I=[[-1]]
   while len(I[0])>0 :

      # 1 if above threshold
      v0   = v[1:-1,1:-1] > threshold

      # Neighbours. 1 if above threshold
      vim1 = numpy.where(v[:-2,1:-1]> threshold,1,0)
      vip1 = numpy.where(v[2:,1:-1] > threshold,1,0)
      vjm1 = numpy.where(v[1:-1,:-2]> threshold,1,0)
      vjp1 = numpy.where(v[1:-1,2:] > threshold,1,0)

      # Sum up neighbours
      tmp = vim1 + vip1 + vjm1 + vjp1

      # If number of neighbours <=num_neighbour, set cell value to threshold
      I = numpy.where(numpy.logical_and(tmp<=num_neighbour, v0  ))
      logger.info("Found %d %d-neighbour cells"%(I[0].size,num_neighbour))
      if I[0].size > 0 :
         #print v[I]
         #print I[0].size,tmp.min(),tmp.max()
         #print v0[I][0], vim1[I][0], vip1[I][0], vjm1[I][0], vjp1[I][0]
         v[1:-1,1:-1][I] = threshold
   return v

def remove_islets(inv,threshold=_default_threshold,num_neighbours=4) :
   """ Remove islets and return modified inv 
   Routine Sets cells with values below threshold to mean of neighbours value if it has more than num_neighbours valid neighbours
   """
   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[0,:] = threshold
   v[:,-1] = threshold
   v[:,0] = threshold
   I=[[-1]]
   while len(I[0])>0 :

      # Points with values below or at threshold
      v0   = v[1:-1,1:-1] <= threshold

      # Neighbouring points above threshold
      fim1 = numpy.where(v[:-2,1:-1]> threshold,1,0)
      fip1 = numpy.where(v[2:,1:-1] > threshold,1,0)
      fjm1 = numpy.where(v[1:-1,:-2]> threshold,1,0)
      fjp1 = numpy.where(v[1:-1,2:] > threshold,1,0)

      # Sum of neighbours above threshold
      tmp = fim1 + fip1 + fjm1 + fjp1
      #I = numpy.where(numpy.logical_and(tmp>=3, v0  ))
      #I = numpy.where(numpy.logical_and(tmp==4, v0  ))
      I = numpy.where(numpy.logical_and(tmp>=num_neighbours, v0  ))
      logger.info("Found %d islets"%I[0].size)

      # Set islet value to mean of neighbours
      if I[0].size > 0 :
         vim1 = v[:-2,1:-1][I] * fim1[I]
         vip1 = v[2:,1:-1] [I] * fip1[I]
         vjm1 = v[1:-1,:-2][I] * fjm1[I]
         vjp1 = v[1:-1,2:] [I] * fjp1[I]
         v[1:-1,1:-1][I] = (vim1 + vip1 + vjm1 + vjp1)/tmp[I]
   return v


def remove_isolated_basins(lon,lat,inv,lon0,lat0,threshold=_default_threshold) :
    """ Find separated regions where inv > threshold.

    If lon0,lat0 is not empty , include regions that contains the points in these arrays
    If lon0,lat0 is     empty , include only the biggest region (in terms of grid cells)
    Returns modified inv, where excluded regions are set to the threshold value
    """
    import scipy.ndimage.measurements

    mask = inv > threshold
    #print numpy.count_nonzero(mask),mask.size
    labarray,num_features=scipy.ndimage.measurements.label(mask)
    logger.info("Found %d features"%num_features)
    feat_count = {}
    for i in range(num_features):
       feat_count[i+1] = numpy.count_nonzero(labarray==i+1)


    # Order by feature count
    tmp = sorted(feat_count.items(), key=lambda x:x[1],reverse=True)
    for i in tmp :
       logger.info( "Feature %03d: %d cells" %(i[0],i[1]))
    main_feature = tmp[0][0] 
    logger.info("Main feature in terms of cells is feature %d"%main_feature)

    # Find points nearest to lon0, lat0
    outv=numpy.ones(inv.shape)*threshold
    if len(lon0)> 0 :
       for lo,la in zip(lon0,lat0) :
         dist = numpy.sqrt((lo-lon)**2+(lat-la)**2) # close enough for this purpose
         I = numpy.argmin(dist)
         if I and inv.flatten()[I] > threshold :
            i,j=numpy.unravel_index(I,inv.shape)
            feature=labarray[i,j]
            logger.info( "Position (%7.3f,%7.3f) : Feature %d is used"%(lo,la,feature))
            outv = numpy.where(labarray==feature,inv,outv)
    else :
       logger.info( "No control points given, using main feature")
       outv = numpy.where(labarray==main_feature,inv,outv)
    return outv

def remove_inconsistent_nesting_zone(depth,threshold=_default_threshold) :
    """ 
    Find points near the boundary where 1st edge cell is land, 2nd edge cell is ocean,  but 3rd and 4th cell out 
    may be land.  python indices start from 0
    """

    # test
    #depth[1,10]=100.
    #depth[2,10]=0.
    #depth[3,10]=0.

    dubious_i2   = depth[:, 1] > threshold
    dubious_iim1 = depth[:,-2] > threshold
    dubious_j2   = depth[ 1,:] > threshold
    dubious_jjm1 = depth[-2,:] > threshold
    for k in range(1,3) :
       dubious_i2     = numpy.logical_and(dubious_i2  ,depth[:,1+k]  <= threshold)
       dubious_iim1   = numpy.logical_and(dubious_iim1,depth[:,-2-k] <= threshold)
       dubious_j2     = numpy.logical_and(dubious_j2  ,depth[1+k,:]  <= threshold)
       dubious_jjm1   = numpy.logical_and(dubious_jjm1,depth[-2-k,:] <= threshold)

    depth[dubious_i2 , 1:4]=threshold
    depth[dubious_iim1,1:4]=threshold
    depth[1:4,dubious_j2  ]=threshold
    depth[1:4,dubious_jjm1]=threshold

    logger.info("Found %d inconsistent nesting points in western  edge"%( numpy.count_nonzero(dubious_i2)))
    logger.info("Found %d inconsistent nesting points in eastern  edge"%( numpy.count_nonzero(dubious_iim1)))
    logger.info("Found %d inconsistent nesting points in southern edge"%( numpy.count_nonzero(dubious_j2)))
    logger.info("Found %d inconsistent nesting points in northern edge"%( numpy.count_nonzero(dubious_jjm1)))

    return depth





def spherdist_haversine(lon1,lat1,lon2,lat2) :
   """ Sphere distance using haversine formula"""
   deg2rad = numpy.pi / 180.
   dlon=lon2-lon1
   dlat=lat2-lat1
   h=hav(dlat*deg2rad) + numpy.cos(lat1*deg2rad) * numpy.cos(lat2*deg2rad) * hav(dlon*deg2rad)
   return 2.* 6371000 * numpy.arcsin(numpy.sqrt(h))
   


def hav(theta) :
   """ Haversine. theta in radians """
   return numpy.sin(theta*.5)**2
   

def p_azimuth(plon1,plat1,plon2,plat2) :
   tmp = fwd_azimuth( plon1,plat1,plon2,plat2)
   # Compass angle -> angle rel lat line
   tmp = numpy.radians(90-tmp)
   # Put in range [-pi,pi)
   tmp = numpy.where(tmp >= numpy.pi,tmp-2*numpy.pi,tmp) # Safest option
   tmp = numpy.where(tmp < -numpy.pi,tmp+2*numpy.pi,tmp) # Safest option
   #mult = max(-tmp.min()/numpy.pi,1)
   #mult = 2*(mult/2) + 1                                # use any odd number as long as tmp+mult*numpy.pi > 0...
   #tmp = numpy.fmod(tmp+9*numpy.pi,2*numpy.pi)-numpy.pi # use any odd number as long as tmp+mult*numpy.pi > 0...
   return tmp

def fwd_azimuth(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return az1


