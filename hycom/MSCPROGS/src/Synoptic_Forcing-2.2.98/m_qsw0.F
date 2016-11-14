      module m_qsw0
      contains
      subroutine qsw0(radfl0,cawdir,clouds,rlat,time,daysinyear)
      use mod_xc
      implicit none

      real, dimension(idm,jdm), intent(out) :: 
     &  cawdir,radfl0
      real, dimension(idm,jdm), intent(in)  :: 
     &  clouds,rlat
      real*8, intent(in) :: time
      integer, intent(in) :: daysinyear

      real pi2,deg,rad,eepsil
      integer ifrac,npart
      real fraci,absh2o,s0,day,dangle,decli,sundv
      real cc,sin2,cos2,cosz,scosz,stot,bioday,biohr,hangle
      real srad,sdir,sdif,altdeg,cfac,ssurf
      integer i,j,l
      real*8 :: daysinyear8

      ! Average number of days in year over a 400-year cycle (Gregorian
      ! Calendar)
      real, parameter :: daysinyear400=365.2425


      daysinyear8=daysinyear


c --- -------------------------------------------------------------------
c --- compute 24 hrs mean solar irrradiance at the marine surface layer
c --- (unit: w/m^2)
c --- -------------------------------------------------------------------

c --- set various quantities

      !print *,'calling qsw0'
      pi2=8.*atan(1.)          !        2 times pi
      deg=360./pi2             !        convert from radians to degrees
      rad=pi2/360.             !        convert from degrees to radians
      eepsil=1.e-9             !        small number

      ifrac=24                 !        split each 12 hrs day into ifrac parts
      fraci=1./ifrac           !        1 over ifrac

      absh2o=0.09              ! ---    absorption of water and ozone
      s0=1365.                 ! w/m^2  solar constant

c --- -------------------------------------------------------------------
c --- compute 24 hrs mean solar radiation at the marine surface layer 
c --- -------------------------------------------------------------------

CKAL  ttime=time+dt/86400.
CKAL? day=aint(time*365./360.)          !accumulated day number (jan1=0,364,..)
      !day=aint(time/float(daysinyear))  !accumulated day number (jan1=0,364,..)

      day=modulo(time,daysinyear8)    !0 < day < 364
      day=floor(day)

      dangle=pi2*day/float(daysinyear)   !day-number-angle, in radians 
      if (day<0. .or. day>daysinyear+1) then
         if (mnproc==1) then
            print *,'qsw0: Error in day for day angle'
            print *,'Day angle is ',day,daysinyear,time
         end if
         call xcstop('(qsw0)')
         stop
      end if

      


c --- compute astronomic quantities -- 
      decli=.006918+.070257*sin(dangle)   -.399912*cos(dangle)
     $             +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)
     $             +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)

      sundv=1.00011+.001280*sin(dangle)   +.034221*cos(dangle)
     $             +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)

c --- scan through each of the grid cells


!$OMP PARALLEL DO PRIVATE(j,l,i,cc,sin2,cos2,scosz,stot,bioday,
!$OMP&     npart,biohr,hangle,srad,sdir,sdif,altdeg,cfac,ssurf,
!$OMP&     cosz)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jdm
      do i=1,idm

c --- compute cloudiness fraction

C KAL cc=clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
C KAL.  +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3
      cc = clouds(i,j)

c --- compute astronomic quantities

      sin2=sin(rlat(i,j))*sin(decli)
      cos2=cos(rlat(i,j))*cos(decli)

c --- split each day into ifrac parts, and compute the solar radiance for 
c --- each part. by assuming symmetry of the irradiance about noon, it
c --- is sufficient to compute the irradiance for the first 12 hrs of
c --- the (24 hrs) day (mean for the first 12 hrs equals then the mean
c --- for the last 12 hrs)

      scosz=0.
      stot=0.

      do npart=1,ifrac
         bioday=day+(npart-.5)*fraci*.5

         biohr=bioday*86400.                !hour of day in seconds
         biohr=amod(biohr+43200.,86400.)    !hour of day;  biohr=0  at noon
         hangle=pi2*biohr/86400.            !hour angle, in radians

         cosz=amax1(0.,sin2+cos2*cos(hangle)) !cosine of the zenith angle
         scosz=scosz+cosz                     !  ..accumulated..

         srad =s0*sundv*cosz                  !extraterrestrial radiation

         !print *,i,j,npart,srad,cosz,eepsil
         ! obs: .7^100 = 3x10^-16 , an already ridicolously low number ...
         sdir=srad*0.7**(min(100.,1./(cosz+eepsil)))    !direct radiation component
!         sdir=srad*0.7**(1./(cosz+eepsil))    !direct radiation component
!         sdir=srad * exp(-0.356674943938732447/(cosz+eepsil))         

         sdif=((1.-absh2o)*srad-sdir)*.5      !diffusive radiation component

         altdeg=amax1(0.,asin(min(1.0,sin2+cos2)))*deg !solar noon altitude in degrees

         cfac=(1.-0.62*cc+0.0019*altdeg)      !cloudiness correction 

         ssurf=(sdir+sdif)*cfac

         stot=stot+ssurf
      enddo

      scosz=scosz*fraci                    !24-hrs mean of  cosz
      radfl0(i,j)=stot*fraci               !24-hrs mean shortw rad in w/m^2

      !cawdir(i,j)=1.-amax1(0.15,0.05/(scosz+0.15)) 
      cawdir(i,j)=1.-amin1(0.15,0.05/(scosz+0.15))  ! Whoa ...
      enddo
      enddo
!$OMP END PARALLEL DO
      !print *,time,day,dangle,decli,sundv,radfl0(63,95)
                                           !co-albedo over water for dir light
      end subroutine qsw0
      end module m_qsw0
