module m_mixlayer_depths

contains

   subroutine mixlayer_depths(tem,sal,dpth,mld1,mld2,idm,jdm,kdm)
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::tem, sal 
      real, dimension(idm,jdm,kdm+1), intent(in)  ::dpth
      real, dimension(idm,jdm)      , intent(out) ::mld1,mld2

      real, parameter :: thref=1e-3
      real :: pi,radian ! Dummys for smt_funcs.h
      include 'stmt_funcs.H'
      integer :: i,j,k,k1,k2

      !real, parameter :: eps_dens=0.05 ;
      real, parameter :: eps_dens=0.03 ;
      real, parameter :: eps_temp=0.2  ;

      real, dimension(kdm) :: tt, ss, dens
      real, dimension(kdm) :: zzmid,dpl
      real, dimension(kdm+1) :: zz
      real tvsz

      do j=1,jdm
      do i=1,idm
      if (dpth(i,j,kdm+1)>10.0) then

         zz=dpth(i,j,1:kdm+1)
         zzmid(:)=(dpth(i,j,2:kdm+1)+dpth(i,j,1:kdm))/2
         tt=tem (i,j,1:kdm)
         ss=sal (i,j,1:kdm)
         do k=1,kdm
            dens(k)=sig(tt(k),ss(k))
            dpl (k)=zz(k+1)-zz(k)
         end do

            
         ! MLD1
         k=1
         do while ( abs(tt(1)-tt(k)) < eps_temp      .and. k<=kdm-1)
            k=k+1
         end do
         k1=k
         !if ( k1 < kdm .and. dpl(k1+1)  > 10.0  .and. &
         !     tt   ( 1) - tt   (k1)  > eps_temp  ) then
         !   tvsz = (tt(k1+1)-tt(k1)) / (zzmid(k1+1) - zzmid(k1)) 
         !   if (tvsz < 1e-10 ) then
         !      mld1(i,j) = (tt(1) - eps_temp - tt(k1)) / tvsz
         !      mld1(i,j) = mld1(i,j) + zzmid(k1)
         !   else
         !      mld1(i,j) =zz(k1)
         !   endif
         !else
         !   mld1(i,j) = zz(k1)
         !end if
         mld1(i,j) = zz(k1)
            

         ! MLD2
         k=1
         do while ( dens(k)-dens(1) < eps_dens  .and. k<=kdm-1)
            k=k+1
         end do
         k2=k
         !if ( k2 < kdm .and. dpl(k2+1)  > 10.0  .and. &
         !     tt   ( 1) - tt   (k2)  > eps_temp  ) then
         !   tvsz = (tt(k2+1)-tt(k1)) / (zzmid(k2+1) - zzmid(k2)) 
         !   if (tvsz < 1e-10 ) then
         !      mld2(i,j) = (tt(1) - eps_temp - tt(k2)) / tvsz
         !      mld2(i,j) = mld1(i,j) + zzmid(k2)
         !   else
         !      mld2(i,j) =zz(k2)
         !   endif
         !else
         !   mld2(i,j) = zz(k2)
         !end if
         mld2(i,j) = zz(k2)
            


         !mld2(i,j) = zz(k)
         !print *,k1,mld1(i,j),k2,mld2(i,j),zz(kdm+1)

      end if
      end do
      end do

   end subroutine
end module








