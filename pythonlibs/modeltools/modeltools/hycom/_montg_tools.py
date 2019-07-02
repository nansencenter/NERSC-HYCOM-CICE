

#     Original fortran code (mod_momtum)
#     # m_prime in lowest layer
#     montg(i,j,kk,1)=psikk(i,j,1)+
#    &        ( p(i,j,kk+1)*(thkk(i,j,1)-thstar(i,j,kk,1))
#    &          -pbavg(i,j,m)*thstar(i,j,kk,1) )*thref**2
#     if     (kapref.eq.-1) then
#             montg(i,j,kk,2)=psikk(i,j,2)+
#    &          ( p(i,j,kk+1)*(thkk(i,j,2)-thstar(i,j,kk,2))
#    &            -pbavg(i,j,m)*thstar(i,j,kk,2) )*thref**2
#           endif !kapref.eq.-1

#     Original fortran code (mod_momtum)
#c ---       m_prime in remaining layers:
#            do k=kk-1,1,-1
#              montg(i,j,k,1)=montg(i,j,k+1,1)+p(i,j,k+1)*oneta(i,j)
#     &            *(thstar(i,j,k+1,1)-thstar(i,j,k,1))*thref**2
#              if     (kapref.eq.-1) then
#                montg(i,j,k,2)=montg(i,j,k+1,2)+p(i,j,k+1)*oneta(i,j)
#     &              *(thstar(i,j,k+1,2)-thstar(i,j,k,2))*thref**2
#              endif !kapref.eq.-1
#            enddo !k

# Example of use:
# ... we have ...
#     montg1 = montg1_no_pb + montg1_pb * pbavg
#     srfhgt = montg1 + thref*pbavg
#            = montg1_no_pb + montg1_pb * pbavg + thref*pbavg
#  ... which gives ...
#     pbavg  = (srfhgt-montg1_no_pb)/(montg1_pb+thref)


import _constants


# Get component of montgomery potential at surface that depends on pbavg
def montg1_pb(thstar,p) :
   kdm=thstar.shape[0]
   thref=_constants.thref
   # Part that depends on pbavg :
   montgpb=-thstar[kdm-1,:,:]*thref**2
   for k in reversed(range(kdm-1)) :
      # PArt that depends on pbavg (through oneta) 
      montgpb[:,:]=montgpb[:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2/p[kdm,:,:]
   return montgpb




# Get component of montgomery potential at surface that does not depend on pbavg
def montg1_no_pb(psikk,thkk,thstar,p) :
   kdm=thstar.shape[0]
   thref=_constants.thref
   montgc=psikk+(p[kdm,:,:]*(thkk-thstar[kdm-1,:,:]))*thref**2
   for k in reversed(range(kdm-1)) :
      montgc [:,:]=montgc [:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2 
   return montgc

