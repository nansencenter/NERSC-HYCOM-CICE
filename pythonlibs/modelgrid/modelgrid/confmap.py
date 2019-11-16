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
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)




# some constants
_pi_1=numpy.pi
_pi_2=_pi_1/2.
_deg=180./_pi_1
_rad=1.0/_deg
_epsil=1.0E-9

#module mod_confmap
#   private
#   real pi_1
#   real pi_2
#   real deg
#   real rad
#   real theta_a
#   real phi_a
#   real theta_b
#   real phi_b
#   real di
#   real dj
#   real dm
#   complex imag
#   complex ac
#   complex bc
#   complex cmna
#   complex cmnb
#   real mu_s
#   real psi_s
#   real epsil
#   logical mercator
#
#   real lat_a,lon_a
#   real lat_b,lon_b
#   real wlim,elim
#   real slim,nlim
#   real mercfac
#   integer ires,jres
#
#   !flags, mainly used by grid generation routine
#   logical :: dotopo, dolatlon, final
#
#   public :: initconfmap, oldtonew, pivotp, bilincoeff, &
#      ll2gind, gind2ll, newtoold, confmap_calctopo,     &
#      getgridP, getgridU, getgridV, getgridQ
#
#contains
#



class ConformalMapping(object) :
   """ Initializes a conformal mapping, described in XXXX"""

   def __init__(self,lat_a,lon_a,lat_b,lon_b,
         wlim,elim,ires,
         slim,nlim,jres,
         mercator,
         mercfac,lold) :
      """Constructor: arguments: 
         lat_a, lon_a     : position of pole A in geo coordinates
         lat_b, lon_b     : position of pole B in geo coordinates
         wlim, elim, ires : western,  eastern  limits in new coords and number of points in 1st dim
         slim, nlim, jres : southern, northern limits in new coords and number of points in 2nd dim
         mercator         : TODO
         mercfac, lold    : TODO
      """


      self._lat_a = lat_a
      self._lon_a = lon_a
      self._lat_b = lat_b
      self._lon_b = lon_b

      self._wlim = wlim
      self._elim = elim
      self._ires = ires

      self._slim = slim
      self._nlim = nlim
      self._jres = jres

      self._mercator = mercator
      self._mercfac  = mercfac 
      self._lold     = lold


      logger.info('Pole A:(%10.4f,%10.4f)'%(self._lat_a,self._lon_a))
      logger.info('Pole B:(%10.4f,%10.4f)'%(self._lat_b,self._lon_b))
      logger.info('W-E  n:(%10.4f,%10.4f,%d)'%(self._wlim,self._elim,self._ires))
      logger.info('S-N  n:(%10.4f,%10.4f,%d)'%(self._slim,self._nlim,self._jres))
      logger.info('MERC  :%r'%self._mercator)
      logger.info('FAC   :%14.6f %r'%(self._mercfac,self._lold))


      self._di=(self._elim-self._wlim)/float(ires-1)   # delta lon'
      self._dj=(self._nlim-self._slim)/float(jres-1)   # delta lat' for spherical grid

      if self._mercator:
         self._dj=self._di
         self._dm=self._di
         if lold:
            logger.info('initconfmap: lold')
            self._slim=-self._mercfac*self._jres*self._dj
         else:
            logger.info('initconfmap: not lold')
            self._slim= self._mercfac

      # transform to spherical coordinates
      self._theta_a=self._lon_a*_rad
      self._phi_a=_pi_2-self._lat_a*_rad
      self._theta_b=self._lon_b*_rad
      self._phi_b=_pi_2-self._lat_b*_rad

      # find the angles of a vector pointing at a point located exactly
      # between the poles
      cx=numpy.cos(self._theta_a)*numpy.sin(self._phi_a)+numpy.cos(self._theta_b)*numpy.sin(self._phi_b)
      cy=numpy.sin(self._theta_a)*numpy.sin(self._phi_a)+numpy.sin(self._theta_b)*numpy.sin(self._phi_b)
      cz=numpy.cos(self._phi_a)+numpy.cos(self._phi_b)

      theta_c=numpy.arctan2(cy,cx)
      self._phi_c=_pi_2-numpy.arctan2(cz,numpy.sqrt(cx*cx+cy*cy))

      # initialize constants used in the conformal mapping
      self._imag=complex(.0,1.)
      self._ac=numpy.tan(.5*self._phi_a)*numpy.exp(self._imag*self._theta_a)
      self._bc=numpy.tan(.5*self._phi_b)*numpy.exp(self._imag*self._theta_b)
      self._c =numpy.tan(.5*self._phi_c)*numpy.exp(self._imag*theta_c)
      self._cmna=self._c-self._ac
      self._cmnb=self._c-self._bc

      w=self._cmnb/self._cmna
      self._mu_s=numpy.arctan2(numpy.imag(w),numpy.real(w))
      self._psi_s=2.*numpy.arctan(abs(w))

        

   def oldtonew(self,lat_o,lon_o) :
      # this routine performes a conformal mapping of the old to the new
      # coordinate system
      lat_o=numpy.atleast_1d(lat_o)
      lon_o=numpy.atleast_1d(lon_o)

      # transform to spherical coordinates
      theta=numpy.mod(lon_o*_rad+3.0*_pi_1,2.0*_pi_1)-_pi_1
      phi=_pi_2-lat_o*_rad

      # transform to the new coordinate system: 1)
      z=numpy.tan(.5*phi)*numpy.exp(self._imag*theta)
      w=(z-self._ac)*self._cmnb/((z-self._bc)*self._cmna)
      mu=numpy.arctan2(numpy.imag(w),numpy.real(w))
      psi=2.*numpy.arctan(abs(w))

      # transform to the new coordinate system: 2)
      I=numpy.abs(phi-_pi_1)<_epsil
      mu[I]=self._mu_s
      psi[I]=self._psi_s

      # transform to the new coordinate system: 3)
      I=numpy.logical_and(numpy.abs(phi-self._phi_b)<_epsil, numpy.abs(theta-self._theta_b)<_epsil)
      mu[I]=0.
      psi[I]=_pi_1

      # transform to new lat/lon coordinates
      lat_n=(_pi_2-psi)*_deg
      lon_n=mu*_deg

      return lat_n,lon_n


   def newtoold(self,lat_n,lon_n) :
      # this routine performes a conformal mapping of the new to the old 
      # coordinate system

      # transform to spherical coordinates
      mu=numpy.mod(lon_n*_rad+3*_pi_1,2*_pi_1)-_pi_1
      psi=numpy.abs(_pi_2-lat_n*_rad)
      
      # transform to the old coordinate system
      w=numpy.tan(.5*psi)*numpy.exp(self._imag*mu)
      z=(self._ac*self._cmnb-w*self._bc*self._cmna)/(self._cmnb-w*self._cmna)
      theta=numpy.arctan2(numpy.imag(z),numpy.real(z))
      phi=2.*numpy.arctan(numpy.abs(z))

      I = numpy.abs(psi-_pi_1) < _epsil
      theta[I]=self._theta_b
      phi  [I]=self._phi_b

      I = numpy.logical_and(
         numpy.abs(mu-self._mu_s)<_epsil,
         (psi-self._psi_s)<_epsil)
      theta[I]=0.
      phi  [I]=_pi_1

      # transform to lat/lon coordinates
      lat_o=(_pi_2-phi)*_deg
      lon_o=theta*_deg

      return lat_o,lon_o


   def pivotp(self,lat_o,lon_o) :
      # This subroutine computes the pivot point of each of the observations
      # in the temporary array tmpobs of type observation. The pivot point
      # is the biggest i and the biggest j, (i,j) is the computation points/
      # the grid, that is less than the position of the observation.
      lat_o=numpy.atleast_1d(lat_o)
      lon_o=numpy.atleast_1d(lon_o)

      # fix for wrap-around
      # western limit in new coordinates can be east of eastern limit (sigh)....
      # in that case di is < 0
      lontmp=numpy.copy(lon_o)
      #
      I=numpy.logical_and(lontmp<self._wlim,self._di>0.)
      lontmp[I]=lontmp[I]+360.
      #
      I=numpy.logical_and(lontmp>self._wlim,self._di<0.)
      lontmp[I]=lontmp[I]+360.
      #
      I=lontmp-self._wlim > 360.
      lontmp[I] = lontmp[I]-360.

      ipiv=(lontmp-self._wlim)/self._di + 1.0

      jpiv=numpy.ones(ipiv.shape)*-999
      if self._mercator: 
         I=numpy.abs(lat_o) < 89.999 
         tmptan=numpy.tan(0.5*_rad*lat_o[I]+0.25*_pi_1)
         jpiv[I]=(numpy.log(tmptan)-self._slim*_rad)/(_rad*self._dj) + 1.0
      else :
         jpiv=(lat_o-self_slim)/dj + 1.0

      # Returns floats, cast to int to use as index/pivot points
      return ipiv,jpiv



   def get_grid_point(self,i,j,shifti=0.,shiftj=0.) :
      #! Used when creating a grid
      #! Retrieves  lon/lat for grid index (i,j) at P points
      #lon_n=self._wlim+(i-1)*self._di+shifti
      #lat_n=self._slim+(j-1)*self._dj+shiftj
      lon_n=self._wlim+(i-1+shifti)*self._di
      lat_n=self._slim+(j-1+shiftj)*self._dj
      if self._mercator :
         lat_n=self._slim+(j-1+shiftj)*self._dm
         lat_n=(2.*numpy.arctan(numpy.exp((lat_n*_rad)))-_pi_1*.5)*_deg
      return self.newtoold(lat_n,lon_n)


   def get_grid(self,shifti=0.0,shiftj=0.,extended=False) :
      #! Used when creating a grid
      #! Retrieves  lon/lat for grid index (i,j) at P points
      if extended :
         i,j=numpy.meshgrid(numpy.arange(self._ires+2),numpy.arange(self._jres+2))
      else  :
         i,j=numpy.meshgrid(numpy.arange(self._ires)+1,numpy.arange(self._jres)+1)
      return self.get_grid_point(i,j,shifti=shifti,shiftj=shiftj)


   #def get_grid_P(self) : return self.get_grid()
   #def get_grid_U(self) : return self.get_grid(shifti=-0.5)
   #def get_grid_V(self) : return self.get_grid(shiftj=-0.5)
   #def get_grid_Q(self) : return self.get_grid(shifti=-0.5,shiftj=-0.5)


   def ll2gind(self,lat_o,lon_o) :
      """ Returns grid index (floats) of specified lat and lon coordinates on the grid"""
      lat_n,lon_n = self.oldtonew(lat_o,lon_o)
      return self.pivotp(lat_n,lon_n)

   def gind2ll(self,i,j,shifti=0.,shiftj=0.):
      """ Returns lat and lon for grid index i,j """
      return self.get_grid_point(i,j,shifti=shifti,shiftj=shiftj)




#subroutine initconfmap(idm,jdm,createin)
#! This routine initialize constants used in the conformal mapping
#! and must be called before the routines 'oldtonew' and 'newtoold'
#! are called. The arguments of this routine are the locations of
#! the two poles in the old coordiante system.
#   implicit none
#   integer, intent(inout) :: idm, jdm
#   logical, intent(in), optional    :: createin


#! local variables
#
#   real cx,cy,cz,theta_c,phi_c
#   complex c,w
#   logical lold,create, ex
#
#   ! If create flag is set, there is no point in comparing with
#   ! idm, jdm. (create may not be the correct name.. it is used by 
#   ! any routine not needing to compare its grid size to the one
#   ! specified by confmap). In that case confmap grid siize is 
#   ! returned in idm, jdm
#   create=.false.
#   if (present(createin)) create=createin

#! Read info file
#   inquire(exist=ex, file='grid.info')
#   if (.not.ex) then
#      print *,'mod_confmap error - grid.info does not exist...'
#      call exit(1)
#   end if
#   open(unit=10,file='grid.info',form='formatted')
#   read(10,*) lat_a,lon_a
#   read(10,*) lat_b,lon_b
#   read(10,*) wlim,elim,ires
#   read(10,*) slim,nlim,jres
#   read(10,*) dotopo
#   read(10,*) dolatlon
#   read(10,*) final
#   read(10,*) mercator
#   read(10,*) mercfac,lold
#   close(10)
#   if (create) then
#      write(*,*)'A:',lat_a,lon_a
#      write(*,*)'B:',lat_b,lon_b
#      write(*,*)'W-E n:',wlim,elim,ires
#      write(*,*)'S-N n:',slim,nlim,jres
#      write(*,*)'TOPO: ',dotopo
#      write(*,*)'LALO: ',dolatlon
#      write(*,*)'FIN: ',final
#      write(*,*)'MERC:',mercator
#      write(*,*)'FAC:', mercfac,lold
#   end if
#   if (.not. create) then
#      if ((ires /= idm).and.(jres /= jdm)) then
#         print *,'initconfmap: WARNING -- the dimensions in grid.info are not'
#         print *,'initconfmap: WARNING -- consistent with idm and jdm'
#         print *,'initconfmap: WARNING -- IGNORE IF RUNNING CURVIINT'
#      endif
#   endif
#
#   ! if create is set, we return grid dimensions 
#   if (create) then
#      idm = ires
#      jdm = jres
#   end if

#   ! some constants
#   pi_1=4.*atan(1.)
#   pi_2=.5*pi_1
#   deg=180./pi_1
#   rad=1.0/deg
#   epsil=1.0E-9
#
#   di=(elim-wlim)/float(ires-1)   ! delta lon'
#   dj=(nlim-slim)/float(jres-1)   ! delta lat' for spherical grid
#
#   if (mercator) then
#      dj=di
#      dm=di
#      if (lold) then
#         print *,'initconfmap: lold'
#         slim=-mercfac*jres*dj
#      else
#         print *,'initconfmap: not lold'
#         slim= mercfac
#      endif
#   endif
#
#! transform to spherical coordinates
#
#   theta_a=lon_a*rad
#   phi_a=pi_2-lat_a*rad
#   theta_b=lon_b*rad
#   phi_b=pi_2-lat_b*rad
#
#! find the angles of a vector pointing at a point located exactly
#! between the poles
#
#   cx=cos(theta_a)*sin(phi_a)+cos(theta_b)*sin(phi_b)
#   cy=sin(theta_a)*sin(phi_a)+sin(theta_b)*sin(phi_b)
#   cz=cos(phi_a)+cos(phi_b)
#
#   theta_c=atan2(cy,cx)
#   phi_c=pi_2-atan2(cz,sqrt(cx*cx+cy*cy))
#
#! initialize constants used in the conformal mapping
#
#   imag=(.0,1.)
#   ac=tan(.5*phi_a)*exp(imag*theta_a)
#   bc=tan(.5*phi_b)*exp(imag*theta_b)
#   c=tan(.5*phi_c)*exp(imag*theta_c)
#   cmna=c-ac
#   cmnb=c-bc
#
#   w=cmnb/cmna
#   mu_s=atan2(aimag(w),real(w))
#   psi_s=2.*atan(abs(w))
#
#end subroutine initconfmap




#subroutine oldtonew(lat_o,lon_o,lat_n,lon_n)
#
#! this routine performes a conformal mapping of the old to the new
#! coordinate system
#
#   implicit none
#
#   real, intent(in)  :: lat_o,lon_o
#   real, intent(out) :: lat_n,lon_n
#
#! local variables
#
#   real theta,phi,psi,mu
#   complex z,w
#
#! transform to spherical coordinates
#
#
#   theta=mod(lon_o*rad+3.0*pi_1,2.0*pi_1)-pi_1
#   phi=pi_2-lat_o*rad
#
#! transform to the new coordinate system
#
#   if (abs(phi-pi_1) < epsil) then
#     mu=mu_s
#     psi=psi_s
#   elseif ((abs(phi-phi_b)<epsil).and.&
#           (abs(theta-theta_b)<epsil)) then
#     mu=.0
#     psi=pi_1
#   else
#     z=tan(.5*phi)*exp(imag*theta)
#     w=(z-ac)*cmnb/((z-bc)*cmna)
#     mu=atan2(aimag(w),real(w))
#     psi=2.*atan(abs(w))
#   endif
#
#! transform to lat/lon coordinates
#
#   lat_n=(pi_2-psi)*deg
#   lon_n=mu*deg
#
#end subroutine oldtonew

#subroutine pivotp(lon,lat,ipiv,jpiv,shiftin)
#! This subroutine computes the pivot point of each of the observations
#! in the temporary array tmpobs of type observation. The pivot point
#! is the biggest i and the biggest j, (i,j) is the computation points/
#! the grid, that is less than the position of the observation.
#
#   implicit none
#
#   real, intent(in) ::  lon,lat
#   integer, intent(out) :: ipiv,jpiv
#   real, intent(in), optional :: shiftin
#
#   real tmptan,tmp
#   real lontmp, shift
#
#   shift=0.
#   if (present(shiftin)) shift=shiftin
#
#! fix for wrap-around
#! western limit in new coordinates can be east of eastern limit (sigh)....
#! in that case di is < 0
#   if (lon < wlim .and. di > 0. ) then
#      lontmp=lon+360.0
#   elseif (lon > wlim .and. di < 0. ) then
#         lontmp=lon-360.0
#!more fixes ...
#   elseif ((lon-wlim)>360.) then
#      lontmp=lon
#      do while(lontmp-wlim>360.)
#         lontmp=lontmp-360.
#      end do
#   else
#      lontmp=lon
#   endif
#
#   ipiv=int((lontmp-wlim)/di)+1
#
#   if (mercator) then
#      if (abs(lat) < 89.999) then
#         tmptan=tan(0.5*rad*lat+0.25*pi_1)
#         !jpiv=int( (log(tmptan)-slim*rad)/(rad*dj) ) +1
#         jpiv=int( (log(tmptan)-slim*rad)/(rad*dj) +1.0 + shift)
#      else
#         jpiv=-999
#      endif 
#   else
#      jpiv=int((lat-slim)/dj + 1.0 + shift)
#   endif
#
#
#! Inverse transformation to check pivot point jpiv
#!   tmp=slim+(jpiv-1)*dj
#!   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg
#
#end subroutine pivotp



#subroutine bilincoeff(glon,glat,nx,ny,lon,lat,ipiv,jpiv,a1,a2,a3,a4)
#! This subroutine uses bilinear interpolation to interpolate the field
#! computed by the model (MICOM) to the position defined by lon, lat
#! The output is the interpolation coeffisients a[1-4]
#! NB  NO locations on land.
#!TODO: periodic grid support
#   implicit none
#
#   integer, intent(in)  :: nx,ny
#   real,    intent(in)  :: glon(nx,ny),glat(nx,ny)
#   real,    intent(in)  :: lon,lat
#   integer, intent(in)  :: ipiv,jpiv
#   real,    intent(out) :: a1,a2,a3,a4
#
#   real t,u
#   real lat_1,lon_1,lat_2,lon_2,lat_n,lon_n,lat_t,lon_t
#
#
#   call oldtonew(glat(ipiv,jpiv),glon(ipiv,jpiv),lat_1,lon_1)
#   call oldtonew(glat(ipiv+1,jpiv+1),glon(ipiv+1,jpiv+1),lat_2,lon_2)
#   call oldtonew(lat,lon,lat_n,lon_n)
#
#!   print *,lat_1,lon_1,lat_2,lon_2,lat_n,lon_n
#!   call oldtonew(glat(ipiv+1,jpiv),glon(ipiv+1,jpiv),lat_t,lon_t)
#!   print *,lat_t,lon_t,lat_t,lon_t
#!
#!   call oldtonew(glat(ipiv,jpiv+1),glon(ipiv,jpiv+1),lat_t,lon_t)
#!   print *,lat_t,lon_t,lat_t,lon_t
#!   print *,'--------------'
#
#
#   t=(lon_n-lon_1)/(lon_2-lon_1)
#   u=(lat_n-lat_1)/(lat_2-lat_1)
#
#   a1=(1-t)*(1-u)
#   a2=t*(1-u)
#   a3=t*u
#   a4=(1-t)*u
#end subroutine bilincoeff
#
#
#

#subroutine gind2ll(ipiv,jpiv,lon,lat)
#! KAL - This routine takes as input floating point ipiv, jpiv 
#! KAL - and calculates actual lon,lat positions from that
#   implicit none
#
#   real, intent(out) ::  lon,lat
#   real, intent(in) :: ipiv,jpiv
#
#   real tmptan,tmp
#   real lontmp
#   real lon_n,lat_n
#
#   !print '(a,2i5,4f10.2)','pivotp:',ipiv,jpiv,lontmp,lon,wlim,di
#   !ipiv=int((lontmp-wlim)/di)+1
#   lon_n = (ipiv-1)*di + wlim
#
#!   if (mercator) then
#!      if (abs(lat) < 89.999) then
#!         tmptan=tan(0.5*rad*lat+0.25*pi_1)
#!         jpiv=int( (log(tmptan)-slim*rad)/(rad*dj) ) +1
#!      else
#!         jpiv=-999
#!      endif 
#!   else
#!      jpiv=int((lat-slim)/dj)+1
#!   endif
#
#   if (mercator) then
#      tmptan = (jpiv-1) * (rad*dj) + slim*rad
#      tmptan = exp(tmptan)
#      lat_n=(atan(tmptan)-0.25*pi_1)*2/rad
#   else
#      lat_n = (jpiv-1) *dj + slim
#   endif
#
#   call newtoold(lat_n,lon_n,lat,lon)
#
#
#! Inverse transformation to check pivot point jpiv
#!   tmp=slim+(jpiv-1)*dj
#!   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg
#
#end subroutine gind2ll

#
#
#
#subroutine ll2gind(lon_o,lat_o,x,y)
#! KAL - This routine takes as input  actual lon, lat
#! KAL - and calculates floating point grid index x,y
#   implicit none
#
#   real, intent(in) ::  lon_o,lat_o
#   real, intent(out) :: x,y
#
#   real tmptan,tmp, lon, lat
#   real lontmp
#
#   call oldtonew(lat_o,lon_o,lat,lon)
#
#! fix for wrap-around
#! western limit in new coordinates can be east of eastern limit (sigh)....
#! in that case di is < 0
#   if (lon < wlim .and. di > 0. ) then
#      lontmp=lon+360.0
#   elseif (lon > wlim .and. di < 0. ) then
#         lontmp=lon-360.0
#!more fixes ...
#   elseif ((lon-wlim)>360.) then
#      lontmp=lon
#      do while(lontmp-wlim>360.)
#         lontmp=lontmp-360.
#      end do
#   else
#      lontmp=lon
#   endif
#
#   !print '(a,2i5,4f10.2)','pivotp:',ipiv,jpiv,lontmp,lon,wlim,di
#   x=(lontmp-wlim)/di+1
#
#   if (mercator) then
#      if (abs(lat) < 89.999) then
#         tmptan=tan(0.5*rad*lat+0.25*pi_1)
#         y=(log(tmptan)-slim*rad)/(rad*dj) +1
#      else
#         y=-999
#      endif 
#   else
#      y=(lat-slim)/dj+1
#   endif
#
#
#! Inverse transformation to check pivot point jpiv
#   tmp=slim+(y-1)*dj
#   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg
#   !print *,'ll2gind',lon_o,lat_o,tmp
#end subroutine ll2gind
#
#
#subroutine newtoold(lat_n,lon_n,lat_o,lon_o)
#! this routine performes a conformal mapping of the new to the old 
#! coordinate system
#
#   implicit none
#
#   real lat_o,lon_o,lat_n,lon_n
#
#! local variables
#
#   real theta,phi,psi,mu
#   complex w,z
#
#! transform to spherical coordinates
#
#   mu=mod(lon_n*rad+3*pi_1,2*pi_1)-pi_1
#   psi=abs(pi_2-lat_n*rad)
#
#! transform to the old coordinate system
#
#   if (abs(psi-pi_1) < epsil) then
#     theta=theta_b
#     phi=phi_b
#   elseif ((abs(mu-mu_s) < epsil).and.((psi-psi_s)<epsil)) then
#     theta=.0
#     phi=pi_1
#   else
#     w=tan(.5*psi)*exp(imag*mu)
#     z=(ac*cmnb-w*bc*cmna)/(cmnb-w*cmna)
#     theta=atan2(aimag(z),real(z))
#     phi=2.*atan(abs(z))
#   endif
#
#! transform to lat/lon coordinates
#
#   lat_o=(pi_2-phi)*deg
#   lon_o=theta*deg
#end subroutine newtoold
#
#
#! Used when creating a grid
#! Retrieves  lon/lat for grid index (i,j) at P points
#subroutine getgridP(i,j,lat_o,lon_o)
#   implicit none
#   integer, intent(in) :: i,j
#   real,    intent(out):: lat_o,lon_o
#   real :: lat_n,lon_n
#   lon_n=wlim+(i-1)*di
#   lat_n=slim+(j-1)*dj
#   if (mercator) then
#      lat_n=slim+(j-1)*dm
#      !print *,i,j,'lat_n',lat_n
#      lat_n=(2.*atan(exp((lat_n*rad)))-pi_1*.5)*deg
#      !print *,i,j,'lat_n',lat_n
#   endif
#   call newtoold(lat_n,lon_n,lat_o,lon_o)
#end subroutine
#
#
#! Used when creating a grid
#! Retrieves  lon/lat for grid index (i,j) at U points
#subroutine getgridU(i,j,lat_o,lon_o)
#   implicit none
#   integer, intent(in) :: i,j
#   real,    intent(out):: lat_o,lon_o
#   real :: lat_n,lon_n
#   ! U point
#   lon_n=wlim+(float(i)-1.0-0.5)*di
#   lat_n=slim+(float(j)-1.0    )*dj
#   if (mercator) then
#      lat_n=slim+(float(j)-1.0   )*dm
#      lat_n=(2.*atan(exp((lat_n*rad)))-pi_1*.5)*deg
#   endif
#   call newtoold(lat_n,lon_n,lat_o,lon_o)
#end subroutine
#
#
#! Used when creating a grid
#! Retrieves  lon/lat for grid index (i,j) at V points
#subroutine getgridV(i,j,lat_o,lon_o)
#   implicit none
#   integer, intent(in) :: i,j
#   real,    intent(out):: lat_o,lon_o
#   real :: lat_n,lon_n
#! V point
#   lon_n=wlim+(float(i)-1.0)*di
#   lat_n=slim+(float(j)-1.0-0.5)*dj
#   if (mercator) then
#      lat_n=slim+(float(j)-1.0-0.5)*dm
#      lat_n=(2.*atan(exp((lat_n*rad)))-pi_1*.5)*deg
#   endif
#   call newtoold(lat_n,lon_n,lat_o,lon_o)
#end subroutine
#
#! Used when creating a grid
#! Retrieves  lon/lat for grid index (i,j) at Q points
#subroutine getgridQ(i,j,lat_o,lon_o)
#   implicit none
#   integer, intent(in) :: i,j
#   real,    intent(out):: lat_o,lon_o
#   real :: lat_n,lon_n
#! q point
#   lon_n=wlim+(float(i)-1.0-0.5)*di
#   lat_n=slim+(float(j)-1.0-0.5)*dj
#   if (mercator) then
#      lat_n=slim+(float(j)-1.0-0.5)*dm
#      lat_n=(2.*atan(exp((lat_n*rad)))-pi_1*.5)*deg
#   endif
#   call newtoold(lat_n,lon_n,lat_o,lon_o)
#end subroutine
#
#! used by grid generation routine
#logical function confmap_calctopo()
#implicit none
#confmap_calctopo = dotopo
#end function
#
#
#
#end module mod_confmap




