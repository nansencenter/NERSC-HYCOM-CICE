      module m_depthave

      contains

      function depthave(field,dpfield,lwdepth,grid)
      use mod_xc
      !use mod_common , ONLY: onem,onemm,depthu,depthv,depths,dp0k
      implicit none
      real,intent(in),dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)::
     &   field, dpfield
      real            , intent(in) :: lwdepth  ! In meters !!
      character(len=1), intent(in) :: grid

      real,dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: depthave

      real,    dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &   dpsum, fdsum, depthx, tp
      integer, dimension(1-nbdy:jdm+nbdy,ms) :: ifx, ilx
      integer, dimension(1-nbdy:jdm+nbdy) :: isx
      real, dimension(itdm,jtdm) :: gdave
      real    :: tmpdp,lwdepth2,updp,lwdepth3
      integer :: i,j,k,l

      include 'common_blocks.h'



      ! Choose grid
      if (grid=='u') then
         isx=isu
         ifx=ifu
         ilx=ilu
         depthx=depthu

      else if (grid=='v') then
         isx=isv
         ifx=ifv
         ilx=ilv
         depthx=depthv

      else if (grid=='p') then
         isx=isp
         ifx=ifp
         ilx=ilp
         depthx=depths

      else
         print *,'Unknown grid id -- '//grid
         stop '(depthave)'
      end if


!$OMP PARALLEL DO PRIVATE(j,l,i,lwdepth2,tmpdp) 
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
           dpsum(i,j)=0.
           fdsum(i,j)=0.
           tp   (i,j)=0.
        enddo
      enddo
!$OMP END PARALLEL DO
      !updp=dp0k(1)
      updp=0.1*onem


      ! Sum up over depth
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,lwdepth2,tmpdp) 
!$OMP&SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isx(j)
         do i=max(1-margin,ifx(j,l)),min(ilx(j,l),ii+margin)

            ! Correct lower depth for bathymetry
            lwdepth2 = min(lwdepth*onem,depthx(i,j)*onem)

            ! Tmpdp is the fraction of this layer included in averaging
            tmpdp    = min(lwdepth2,tp(i,j)+dpfield(i,j,k)) 
     &         - max(updp,tp(i,j)) 
            tmpdp    = max(0.,tmpdp)
            !print *,i,j,k,tmpdp/onem

            ! Cumulative depth tp, dpsum and fdsum is cumulative down to 
            ! max depth (lwdepth)
            tp(i,j)=tp(i,j)+dpfield(i,j,k)
            dpsum(i,j) = dpsum(i,j) + tmpdp
            fdsum(i,j) = fdsum(i,j) + tmpdp*field(i,j,k)
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
cdiag print *,maxval(dpsum)/onem,minval(dpsum,mask=dpsum>onem)/onem
cdiag print *,maxval(fdsum)/(lwdepth*onem),minval(fdsum)/(lwdepth*onem)
cdiag print *,maxval(field)



      ! Average over depth
      depthave=0.
!$OMP PARALLEL DO PRIVATE(j,l,i) 
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isx(j)
      do i=max(1-margin,ifx(j,l)),min(ilx(j,l),ii+margin)
         depthave(i,j) = fdsum(i,j) / (dpsum(i,j) + onemm)
      end do
      end do
      end do
!$OMP END PARALLEL DO
cdiag print *,maxval(depthave)
cdiag print *


      !Test
      !call xcaget(gdave,depthave,0)
      !open(10,file='depthave'//grid//'.tec')
      !WRITE(10,*)'TITLE=""'
      !WRITE(10,*)'VARIABLES=i,j,fld'
      !WRITE(10,*)'ZONE I=',itdm,',J=',jtdm,',F=BLOCK'
      !WRITE(10,'(30i4)')((i,i=1,itdm),j=1,jtdm)
      !WRITE(10,'(30i4)')((j,i=1,itdm),j=1,jtdm)
      !WRITE(10,'(10e15.6)')((gdave(i,j),i=1,itdm),j=1,jtdm)
      !close(10)


      end function depthave
      end module m_depthave





