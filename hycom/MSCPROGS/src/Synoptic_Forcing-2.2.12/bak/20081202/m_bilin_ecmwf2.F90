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

   integer i,j,ia,ib,ja,jb,ifalse,maxlati,minlati,jposb
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid
   integer :: numerr, allerr
   real    :: newlon2

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

   lerr=.false.
   numerr=0
   ! Start interpolation
   !$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jposb,lon,ib,aa,bb,a1,a2,a3,a4,newlon2) &
   !$OMP SHARED (lerr) SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
   do j=1,jdm
   do i=1,idm

         newlon2 = mod(newlon(i,j)+360.d0,360.d0)

         ipos=int((newlon2     -olonref)/odlon+1.0)
         jpos=int((newlat (i,j)-olatref)/odlat+1.0)
         if (jpos>ony) then
            write(*,'(3i5,3f10.2)')jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
            print *,'bilin_ecmwf2: jpos',jpos
            stop '(bilin_ecmwf2)'
         end if
         jposb=min(jpos+1,ony)
         lon=olonref+(ipos-1)*odlon

         ! Assumes periodic data set
         ib=mod(ipos,onx)+1

         aa=(newlon2      - olonref-float(ipos-1)*odlon)/odlon
         bb=(newlat (i,j) - olatref-float(jpos-1)*odlat)/odlat

         if ((aa > 1.0).or.(aa < 0.0)) then
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2     ,lon+odlon
            print *,'bilin_ecmwf2: invalid aa',aa
            stop '(bilin_ecmwf2)'
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            write(*,'(3i5,3f10.2)')i,j,jpos,odlat*(jpos-1),newlat(i,j),olatref+odlat*jpos
            print *,'bilin_ecmwf2: invalid bb',bb
            bb=max(0.0,min(0.999,bb))
            stop '(bilin_ecmwf2)'
         endif

         ! Error in case of wrong aa or bb
         numerr = numerr + abs(floor(aa))
         numerr = numerr + abs(floor(bb))

         a1=(1.0-aa)*(1.0-bb)
         a2=aa*(1.0-bb)
         a3=aa*bb
         a4=(1.0-aa)*bb

         new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jposb)+a4*old(ipos,jposb)
      !endif
   enddo
   enddo
   !$OMP END PARALLEL DO


end subroutine 
end module m_bilin_ecmwf2
