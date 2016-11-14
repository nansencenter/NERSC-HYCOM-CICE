module m_bilin_ecmwf2
contains
subroutine bilin_ecmwf2(old,onx,ony,olonref,olatref,odlon,odlat,new,newlon,newlat,depths)   
   use mod_xc
   implicit none
   integer, intent(in) :: onx         ! x-dimension of old field
   integer, intent(in) :: ony         ! y-dimension of old field
   real, intent(in)    :: old(onx,ony)! old grid
   real, intent(in)    :: olonref     ! Lon - reference point
   real, intent(in)    :: olatref     ! Lat - reference point
   real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
   real, intent(in)    :: odlat       ! Lat - grid spacing in old grid

   real, intent(out)   :: new   (idm,jdm)  ! New interpolated field
   real, intent(in)    :: newlon(idm,jdm) ! Longitudes for new grid
   real, intent(in)    :: newlat(idm,jdm) ! Latitudes for new grid
   real, intent(in)    :: depths(idm,jdm) ! Latitudes for new grid




   integer i,j,ia,ib,ja,jb,ifalse,maxlati,minlati
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid
   integer :: numerr, allerr
   real    :: newlon2(idm,jdm)

   real aa,bb,a1,a2,a3,a4,maxlat,lon,minlat
   character(len=2) tag2
   logical :: lerr




   maxlat=olatref+(ony-1)*odlat
   maxlat=max(olatref,maxlat) ; ! For odlat<0
   minlat=olatref+(ony-1)*odlat
   minlat=min(minlat,olatref)   ! For odlat<0

   ! Set "minimumlat" and "maximumlat" index..
   if (olatref<olatref+(ony-1)*odlat) then
      minlati=1
      maxlati=ony
   else
      minlati=ony
      maxlati=1
   end if

   ! Correct newlon to be in range of  "old" data
   newlon2=newlon
   !print *,'bilin_ecmwf2: max min lon:',maxval(newlon2),minval(newlon2)
   if (any(newlon2<olonref).or.any(newlon2>olonref+360.)) then
      !print *,'Warning - bilin_ecmwf2 : Correcting longitudes'
      where (newlon2<olonref     ) newlon2 = newlon2 + 360.
      where (newlon2>olonref+360.) newlon2 = newlon2 - 360.
   end if
      



   lerr=.false.
   numerr=0
   ! Start interpolation
   !$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,lon,ib,aa,bb,a1,a2,a3,a4) &
   !$OMP SHARED (lerr) SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
   do j=1,jdm
   do i=1,idm

         ipos=int((newlon2(i,j)-olonref)/odlon+1.0)
         jpos=int((newlat(i,j)-olatref)/odlat+1.0)
!        jpos=min(max(1,jpos),ony-1)
         jpos=min(max(1,jpos),ony)
         jb  =min(max(1,jpos+1),ony)

         ! Error in case out of data domain
         !numerr = numerr + (1-sign(1,jpos-1))/2
         !numerr = numerr + (1-sign(1,ipos-1))/2
         !numerr = numerr + (1-sign(1,onx-ipos))/2
         !numerr = numerr + (1-sign(1,ony-(jpos-1)))/2

         lon=olonref+(ipos-1)*odlon

         ! Assumes periodic data set
         ib=mod(ipos,onx)+1

         aa=(newlon2(i,j) - olonref-float(ipos-1)*odlon)/odlon
         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat

         if ((aa > 1.0).or.(aa < 0.0)) then
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2(i,j),lon+odlon
            print *,'bilin_ecmwf2: invalid aa',aa
            lerr=.true.
            numerr=numerr+1
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            print *,'bilin_ecmwf2: invalid bb',bb
            bb=max(0.0,min(0.999,bb))
            lerr=.true.
            numerr=numerr+1
         endif

         a1=(1.0-aa)*(1.0-bb)
         a2=aa*(1.0-bb)
         a3=aa*bb
         a4=(1.0-aa)*bb

         new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jb)+a4*old(ipos,jb)
      !endif
   enddo
   enddo
   !$OMP END PARALLEL DO

   allerr=numerr
   

   if (numerr>0) then
      write(lp,*) 'An error occured in bilin_ecmwf2..'
      write(lp,'(a,i6,a)') 'An error occured in bilin_ecmwf2..', allerr,' errors'
      call flush(lp)
      call exit(1)
   end if


end subroutine 
end module m_bilin_ecmwf2
