import logging
import numpy

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



def isopycnal_coordinate_layers(dprofi,sigprof,mindp,sigma,isotop) :
   # dprof    [nz+1]   : depth profile in orignial data. sigprof of layer i occupies depth range(dprof[i],dprof[i+1])
   # sigprof  [nz,np]  : density profile in orignial data
   # mindp    [np,kdm] : min layer thickness in original data
   #
   #TODO: Make sure vertical levels rank match
   nz=sigprof.shape[0]
   kdm=sigma.size

   # MAx depth in data. TODO. Get from profiles
   maxd=numpy.zeros(sigprof[0,:].shape)
   for k in range(nz) :
      maxd[~sigprof[k,:].mask] = dprofi[k]
   logger.info("MAx depth in data is %d m"%numpy.max(maxd))

   # Loop over output layers
   newintf=numpy.zeros((mindp.shape[0],kdm+1))  # New layer interfaces
   intsig=numpy.zeros(mindp.shape)              # Integrated sigma value
   for k in range(kdm) :

      # Target layer upper interface
      upint=newintf[:,k]

      # Mix water over integration range
      dpsum=numpy.zeros(upint.shape)
      sgsum=numpy.zeros(upint.shape)
      sg   =numpy.zeros(upint.shape)
      for k2 in range(nz) :


         # Range of this layer
         upint2=dprofi[k2]
         lwint2=dprofi[k2+1]

         # Part of this layer in target layer
         u = numpy.maximum(upint2,upint)
         l = lwint2
         dp=numpy.maximum(0.,l-u)

         #Mask where dpsum > 0  
         Imask = dpsum > 0.

         ####################  dpsum > 0  and summed layer density < target density ###########

         # Integrated value of sigma up to this point
         sg[Imask]=sgsum[Imask]/dpsum[Imask]

         # Mask where target layer heavier than integrated value => Ok to add more layers
         Jmask = numpy.logical_and(Imask,sigma[k] > sg)
         #Jmask = numpy.logical_and(Jmask,~salprof.mask[k2,:])
         Jmask = numpy.logical_and(Jmask,~sigprof.mask[k2,:])

         # Fraction of layer to be added
         dpfrac = (dpsum[Jmask]*sg[Jmask] - dpsum[Jmask]*sigma[k]) / (sigma[k] - sigprof[k2,Jmask])

         # Can not add more than dp!
         dpfrac = numpy.minimum(dpfrac,dp[Jmask]) #

         # Dpfrac can be negative if sigprof[k2,Jmask] < sigma[k] ( must mix "negatively"). 
         # In that case add entire layer. Perhaps next heavy layer can tip the scale ...
         dpfrac=numpy.where(dpfrac < 0.,dp[Jmask],dpfrac)

         # Update values 
         dpsum[Jmask]=dpsum[Jmask]+dpfrac
         sgsum[Jmask]=sgsum[Jmask]+sigprof[k2,Jmask]*dpfrac

         # NB: No need to treat summed layer density > target density, since mixing wont 
         # enable us to reach target density. (sg can only increase as we sum deeper)


         ####################  dpsum == 0  and sigma[k] > sigprof[k2,:] ######################

         #sg = sigprof[k2,k2]

         # target layer heavier than layer value. Ok to add more layers to mix downward
         Jmask = numpy.logical_and(~Imask,sigma[k] >  sigprof[k2,:])
         #Jmask = numpy.logical_and(Jmask,~salprof.mask[k2,:])
         Jmask = numpy.logical_and(Jmask,~sigprof.mask[k2,:])
         dpsum[Jmask]=dpsum[Jmask]+dp[Jmask]
         sgsum[Jmask]=sgsum[Jmask]+sigprof[k2,Jmask]*dp[Jmask]

         # No need to treat new layer density > target layer density, since mixing wont enable 
         # us to reach target density

         #################### ######################################### ######################




      # end k2 loop
      #print dpsum

      # If layer below isotop, we must have fixed coords
      dpsum=numpy.where(newintf[:,k] < isotop,mindp[:,k],dpsum)

      # Make sure dpsum adheres to minimum layer thickness
      dpsum=numpy.maximum(dpsum,mindp[:,k])
      
      # Adjust layer interface
      newintf[:,k+1] = newintf[:,k] + dpsum

      # MAke sure lowest interface is above sea floor
      newintf[:,k+1] = numpy.minimum(newintf[:,k+1],maxd)
      
      # Effective layer thickness
      dpsum = newintf[:,k+1] - newintf[:,k]

      # Estimated sigma
      msk=dpsum>0.
      intsig[msk,k]=sgsum[msk]/dpsum[msk]

      # TODO: Keep track of final layer density for all layers (not just isopycnals)
   # End k loop


   return newintf,intsig
