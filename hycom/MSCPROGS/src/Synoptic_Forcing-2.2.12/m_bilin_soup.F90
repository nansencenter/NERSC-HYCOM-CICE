module m_bilin_soup
contains
subroutine bilin_soup(old,onx,ony,olonref,olatref,odlon,odlat,new,nx,ny,newlon,newlat)   
   use mod_xc, only :  xcstop, mnproc
   implicit none
   integer, intent(in) :: onx         ! x-dimension of old field
   integer, intent(in) :: ony         ! y-dimension of old field
   real, intent(in)    :: old(onx,ony)! old grid
   real, intent(in)    :: olonref     ! Lon - reference point
   real, intent(in)    :: olatref     ! Lat - reference point
   real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
   real, intent(in)    :: odlat       ! Lat - grid spacing in old grid

   integer, intent(in) :: nx          ! x-dimension of old field
   integer, intent(in) :: ny          ! y-dimension of old field
   real, intent(out)   :: new(nx,ny)  ! New interpolated field
   real, intent(in)    :: newlon(nx,ny) ! Longitudes for new grid
   real, intent(in)    :: newlat(nx,ny) ! Latitudes for new grid




   integer i,j,ia,ib,ja,jb,ifalse
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid
   integer countgrid

   real aa,bb,a1,a2,a3,a4,maxlat,lon
   character(len=2) tag2
   logical :: lerr

!   print *,'XXX',olonref,odlon,onx
!   print *,'XXX',olatref,odlat,ony


   maxlat=olatref+(ony-1)*odlat
   lerr=.false.
   countgrid=0
   ! Start interpolation
   !$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,lon,ib,aa,bb,a1,a2,a3,a4) &
   do j=1,ny
   do i=1,nx
      if (newlat(i,j) < olatref) then
         !KALprint *,'lat < latref'
         new(i,j)=0.0
      elseif (newlat(i,j) > maxlat) then
         ipos=int((newlon(i,j)-olonref)/odlon+1.0)
         new(i,j)=0.0
         !KALprint '(a,2i5,2f8.2,g10.2)','lat > maxlat',i,j,newlat(i,j),newlon(i,j),new(i,j)
      else
         ipos=int((newlon(i,j)-olonref)/odlon+1.0)
         jpos=int((newlat(i,j)-olatref)/odlat+1.0)
         if ((ipos < 1).or.(ipos > onx)) then
            !print *,'bilin_ecmwf: ipos=',ipos,onx,newlon(i,j),olonref
            !stop
            !lerr=.true.
            new(i,j)=0.0
            cycle
         endif
         if ((jpos < 1).or.(jpos > ony)) then
            !print *,'bilin_ecmwf: jpos=',jpos,ony,newlat(i,j),olatref
            !lerr=.true.
            !stop
            new(i,j)=0.0
            cycle
         endif

         countgrid=countgrid+1

         lon=olonref+(ipos-1)*odlon

         if (ipos==onx) then
               ib=1
         else
               ib=ipos+1
         endif

         aa=(newlon(i,j) - olonref-float(ipos-1)*odlon)/odlon
         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat

         if ((aa > 1.0).or.(aa < 0.0)) then
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
            print *,'bilin_ecmwf: invalid aa',aa
            !stop
            lerr=.true.
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            print *,'bilin_ecmwf: invalid bb',bb
            lerr=.true.
            !stop
         endif

         a1=(1.0-aa)*(1.0-bb)
         a2=aa*(1.0-bb)
         a3=aa*bb
         a4=(1.0-aa)*bb

         new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jpos+1)+a4*old(ipos,jpos+1)
      endif
   enddo
   enddo
   !$OMP END PARALLEL DO

   if (lerr) then
      write(*,*) 'An error occured in bilin_ecmwf..'
      call xcstop('(m_bilin_ecmwf.F90)')
      stop '(m_bilin_ecmwf.F90)'
   end if
   
   if (mnproc==1) then
      write(*,*) 'bilin_soup: ',countgrid , ' model points on soup grid',mnproc,'prc. nr.'
   end if

end subroutine 
end module m_bilin_soup
