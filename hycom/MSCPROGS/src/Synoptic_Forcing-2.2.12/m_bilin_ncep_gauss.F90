module m_bilin_ncep_gauss
contains
subroutine bilin_ncep_gauss(old,onx,ony,olon,olat,new,newlon,newlat)   
   use mod_xc
   implicit none
   integer, intent(in) :: onx         ! x-dimension of old field
   integer, intent(in) :: ony         ! y-dimension of old field
   real, intent(in)    :: old(onx,ony)! old grid
   real, intent(in)    :: olon(onx)! old grid
   real, intent(in)    :: olat(ony)! old grid

   real, intent(out)   :: new   (idm,jdm)  ! New interpolated field
   real, intent(in)    :: newlon(idm,jdm) ! Longitudes for new grid
   real, intent(in)    :: newlat(idm,jdm) ! Latitudes for new grid

   integer i,j,ia,ib,ja,jb,ifalse,io
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid

   real aa,bb,a1,a2,a3,a4,lon,olonref,odlon,odlat
   character(len=2) tag2
   real :: valvec(4),mskvec(4)
   integer :: nsea

   olonref= minval(olon)
   odlon   = (olon(2)-olon(1))

! Start interpolation
   do j=1,jdm
   do i=1,idm

      ipos=int((newlon(i,j)-olonref)/odlon+1.0)
      ib=mod(ipos,onx)+1

      ! Substitute formula for "gaussian" grid here 
      if (olat(ony)<olat(1)) then

         jpos=1
         do io=ony,1,-1
            if (olat(io)<newlat(i,j)) then
               jpos=io
            end if
         end do
         jb=max(jpos-1,1)
         if (jpos==1) then
            odlat=(90.-olat(jpos))
         else 
            odlat=(olat(jpos-1)-olat(jpos))
         end if
      else
         print *,'Gaussian grid for increasing lat needs testing'
         stop
         jpos=ony
         do io=1,ony
            if (olat(io)<newlat(i,j)) then
               jpos=io
            end if
         end do
         jpos=min(max(1,jpos),ony)
         odlat=(olat(jpos-1)-olat(jpos))
      end if


      aa=(newlon(i,j) - olon(ipos))/odlon
      bb=(newlat(i,j) - olat(jpos))/odlat

      ! Error in case of wrong aa or bb
      if (aa>1.0 .or. aa<0.0) then
         print *,'Error in bilin_ncep_gauss, aa=',aa
         stop '(bilin_ncep_gauss)'
      elseif (bb>1.0 .or. bb<0.0) then
         print *,'Error in bilin_ncep_gauss, bb=',bb
         print *,newlat(i,j),ipos,jpos,olat(jpos), olat(jb)
         stop '(bilin_ncep_gauss)'
      end if

      a1=(1.0-aa)*(1.0-bb)
      a2=aa*(1.0-bb)
      a3=aa*bb
      a4=(1.0-aa)*bb

      new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jb)+a4*old(ipos,jb)

   enddo
   enddo

end subroutine 
end module m_bilin_ncep_gauss
