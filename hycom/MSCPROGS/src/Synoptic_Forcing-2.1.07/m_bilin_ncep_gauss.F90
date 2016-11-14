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






   integer i,j,ia,ib,ja,jb,ifalse,io,numerr,allerr
   integer ipos    !  index of i-pivot grid point in old grid
   integer jpos    !  index of j-pivot grid point in old grid

   real aa,bb,a1,a2,a3,a4,maxlat,lon,minlat,olonref,odlon,odlat
   character(len=2) tag2
   real :: valvec(4),mskvec(4)
   integer :: nsea


   maxlat=maxval(olat)
   minlat=minval(olat)
   olonref= minval(olon)
   odlon   = (olon(2)-olon(1))
   numerr=0

! Start interpolation
   do j=1,jdm
   do i=1,idm

      ipos=int((newlon(i,j)-olonref)/odlon+1.0)

      ! Substitute formula for "gaussian" grid here 
      jpos=1
      do io=ony,1,-1
         if (olat(io)<newlat(i,j)) then
            jpos=io
         else
            exit
         end if
      end do
      jpos=min(max(2,jpos),ony)
      odlat=(olat(jpos-1)-olat(jpos))

      ! Error in case out of data domain
      numerr = numerr + (1-sign(jpos-1,1))/2
      numerr = numerr + (1-sign(ipos-1,1))/2
      numerr = numerr + (1-sign(onx-ipos,1))/2
      numerr = numerr + (1-sign(ony-jpos,1))/2

      ! Wrap-around if ipos==nx
      ib=mod(ipos,onx)+1
      jb=min(jpos+1,ony)

      aa=(newlon(i,j) - olon(ipos))/odlon
      bb=(newlat(i,j) - olat(jpos))/odlat

      ! Error in case of wrong aa or bb
      numerr = numerr + abs(floor(aa))
      numerr = numerr + abs(floor(bb))

      a1=(1.0-aa)*(1.0-bb)
      a2=aa*(1.0-bb)
      a3=aa*bb
      a4=(1.0-aa)*bb

      new(i,j) = a1*old(ipos,jpos)+a2*old(ib,jpos)+a3*old(ib,jpos-1)+a4*old(ipos,jpos-1)

   enddo
   enddo


   !print *,newlon(1,:)
   

   if (allerr>0) then
      if (mnproc==1) then
         write(lp,*) 'An error occured in bilin_ncep_gauss..'
         write(lp,'(a,i6,a)') 'An error occured in bilin_ncep_gauss..', allerr,'errors'
         call flush(lp)
      end if
      call xcstop('(m_bilin_ncep_gauss.F90)')
      stop '(m_bilin_ncep_gauss.F90)'
   end if

end subroutine 
end module m_bilin_ncep_gauss
