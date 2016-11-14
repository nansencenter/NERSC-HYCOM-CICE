program topography
   implicit none

   integer, parameter :: nx=181     ! i-dim of model grid
   integer, parameter :: ny=240     ! j-dim of model grid
   real depths(1:nx,1:ny)           ! model depths
   real modlon(1:nx,1:ny)           ! model longitudes
   real modlat(1:nx,1:ny)           ! model latitudes

   integer, parameter :: nrx=4320   ! xdim of ETOPO5
   integer, parameter :: nry=2160   ! ydim of ETOPO5
   integer ivecf(nrx,nry)           ! ETOPO5 file


   integer i,j,ix,jy
   real s1,s2,s3,s4,aa,bb
   real, parameter :: five_min=1.0/12.0  
      
   character(len=7)tag7

   integer shorder
   integer nrfilt
   integer shdim
   parameter(shdim=16)
   real sh(0:shdim)



   logical mask(1:nx,1:ny)
   real depths2(1:nx,1:ny)
   real depths3(1:nx,1:ny)


! Reading lat-lon for model gridpoints.
   open(10,file='newpos.dat')
      read(10,'(10g17.10)')modlat,modlon
   close(10)
   write(*,*)'Model lat-lon is read'

! Reading ETOPO5 data set
   open(10,file='ETOPO5.uf')
      read(10)ivecf
   close(10)
   write(*,*)'ETOPO5 is read'

        
! Bilinear interpolation
   do j=1,ny
      do i=1,nx
         ix=int((modlon(i,j))*12.0+1.0)
         jy=int((90.0-modlat(i,j))*12.0+2.0)

         depths(i,j)=float(ivecf(ix,jy))

         ix=max(ix,1)
         ix=min(ix,nrx)
         jy=max(jy,1)
         jy=min(jy,nry)

         s1=float(ivecf(ix,jy)    )
         s2=float(ivecf(ix+1,jy)  )
         s3=float(ivecf(ix+1,jy-1))
         s4=float(ivecf(ix,jy-1)  )

         aa=( modlon(i,j)-float(ix-1)*five_min  )/five_min
         bb=( modlat(i,j)-(90.0-float(jy-1)*five_min) )/five_min
!         write(21,'(2i4,4f10.4)')i,j,float(ix-1)*five_min ,  modlon(i,j), float(ix)*five_min,aa
!         write(21,'(2i4,4f10.4)')i,j,90.0-float(jy-1)*five_min,  modlat(i,j), 90.0-float(jy-2)*five_min,bb
!         write(21,*)
         depths(i,j)=(1.0-aa)*(1.0-bb)*s1+aa*(1.0-bb)*s2+aa*bb*s3+(1.0-aa)*bb*s4
      enddo
   enddo

   mask=.true.
   where (depths >= -15.0) 
      depths=0.0
      mask=.false.
   endwhere

   depths=-depths

! ANITA...
   depths(1:32,:)=0.0        !northern polar regions
   depths(36:53,6:28)=0.0    !Baltic 
   depths(30:37,14:28)=0.0    !Russian lakes
   depths(55:76,1:39)=0.0    !mediterranian, red sea,caspian sea
   depths(56:78,36:40)=0.0   !Persian gulf
   depths(76:85,24:34)=0.0   !Persian gulf
   depths(31:51,177:193)=0.0 !Hudson Bay
   depths(53:56,196:200)=0.0 !    " area
   depths(52:55,200:202)=0.0 !    "
   depths(64:66,237:240)=0.0 !    "
   depths(47:54,237:240)=0.0 !    "
   depths(45:49,236:238)=0.0 !Mediterranian, North atlantic
   depths(44:44,238:238)=0.0 !    "
!   depths(55:59,238:240)=0.0 !    "
!   depths(51:67,86:93)=0.0   !Sea of japan
!   depths(51:66,93:95)=0.0   !    "
!   depths(56:58,95:96)=0.0   !    "
!   depths(40:53,104:104)=0.0 !sea of ohktosk
!   depths(39:50,105:107)=0.0 !    "
!   depths(40:48,107:109)=0.0 !    "
!   depths(38:44,108:110)=0.0 !    "
!   depths(31:37,126:132)=0.0 !Berind strait
   depths(65:66,162:163)=0.0 !Gulf of california
   depths(69:73,164:167)=0.0 !      "
   depths(73:75,166:168)=0.0 !      "

!   shorder=8
   depths2=depths
!   do i=1,5
!      call shfact(shorder,sh,shdim)
!      call shfilt2(shorder,sh,depths2,depths3,nx,ny,shdim)
!   enddo
!   where (.not.mask) depths2=0.0

   write(tag7,'(i3.3,a,i3.3)')nx,'x',ny
   call tecfld('tmp'//tag7,nx,ny,modlat,modlon,depths,depths2)
   open(10,file='depths'//tag7//'.dat')
      write(10,'(10g17.10)')depths
   close(10)

end program

subroutine tecfld(fname,ii,jj,modlat,modlon,depths,field)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj
   real, intent(in) :: modlat(ii,jj)
   real, intent(in) :: modlon(ii,jj)
   real, intent(in) :: depths(ii,jj)
   real, intent(in) :: field(ii,jj)

   character*4 tag

   open(44,file=fname//'.dat',status='unknown')
      write(44,*)'TITLE = "',fname,'"'
      write(44,*)'VARIABLES = "i"  "j"  "lat" "lon" "depths" "depths2"'
      write(44,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in Grid index
      write(44,'(30I4)')((i,i=1,ii),j=1,jj)
      write(44,'(30I4)')((j,i=1,ii),j=1,jj)


! Positions in Lon Lat
      write(44,900)((modlat(i,j),i=1,ii),j=1,jj)
      write(44,900)((modlon(i,j),i=1,ii),j=1,jj)
      write(44,900)((depths(i,j),i=1,ii),j=1,jj)
      write(44,900)((field(i,j),i=1,ii),j=1,jj)

   close(44)
900  format(10(1x,e12.5))

end subroutine tecfld



      subroutine shfact(n,sh,shdim)
      integer n
      integer i,j,shdim
      real sh(0:shdim)
      real f2n,fj,f2nj,ff
      if(n.GT.8)then
      write(*,*)'Error in "shfact"'
      write(*,*)'n is to large (>8): n=',n
      stop
      endif
      if((n.EQ.1).OR.(n.EQ.2).OR.(n.EQ.4).OR.(n.EQ.8))then
      ff=2.0**(2*n)
      f2n=1
      do j=1,2*n
      f2n=f2n*float(j)
      enddo
      fj=1
      sh(0)=((-1.0)**(n-1))/ff
      do j=1,n
      
      f2nj=1
      do i=1,2*n-j
      f2nj=f2nj*float(i)
      enddo
      
      fj=fj*float(j)
      
      sh(j)=(-1.0)**(n+j-1)*f2n/(ff*fj*f2nj)
      enddo
      else
      write(*,*)'error in shfact.  n=',n
      stop
      endif
      return
      end
      subroutine shfilt2(ish,sh,x,y,nx,ny,shdim)
      integer j,i,m,mm
      integer ish,shdim
      real sh(0:shdim)
      real x(nx,ny)
      real y(nx,ny)
      if((2*ish+1.GT.nx).OR.(2*ish+1.GT.ny))then
      write(*,*)'shfilt2:  The domain is to small for ish=',ish
      return
      endif
      
      
      call scopy(nx*ny,0.0,0,y,1)
      
      do i=1,nx
      y(i,1)=x(i,1)
      y(i,ny)=x(i,ny)
      end do
      
      do m=2,ish
      do j=0,ish-1
      mm=1-(m-ish+j)
      if(mm.GE.0)then
      do i=1,nx
      y(i,m)=y(i,m)+sh(j)*(2.0*x(i,1)-x(i,1+mm)+x(i,m+ish-j))
      end do
      else
      do i=1,nx
      y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
      end do
      endif
      enddo
      do i=1,nx
      y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
      end do
      enddo
      
      do m=ish+1,ny-ish
      do j=0,ish-1
      do i=1,nx
      y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
      end do
      enddo
      do i=1,nx
      y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
      end do
      enddo
      
      do m=ny-ish+1,ny-1
      do j=0,ish-1
      mm=(m+ish-j)-ny
      if(mm.GE.0)then
      do i=1,nx
      y(i,m)=y(i,m)+sh(j)*(2.0*x(i,ny)-x(i,ny-mm)+x(i,m-ish+j))
      end do
      else
      do i=1,nx
      y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
      end do
      endif
      enddo
      do i=1,nx
      y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
      end do
      enddo
      
      
      call scopy(nx*ny,0.0,0,x,1)
      
      do m=1,ny
      x(1,m)=y(1,m)
      x(nx,m)=y(nx,m)
      end do
      
      do i=2,ish
      do j=0,ish-1
      mm=1-(i-ish+j)
      if(mm.GE.0)then
      do m=1,ny
      x(i,m)=x(i,m)+sh(j)*(2.0*y(1,m)-y(1+mm,m)+y(i+ish-j,m))
      end do
      else
      do m=1,ny
      x(i,m)=x(i,m)+sh(j)*(y(i+ish-j,m)+y(i-ish+j,m))
      end do
      endif
      enddo
      do m=1,ny
      x(i,m)=(1.0+sh(ish))*y(i,m)+x(i,m)
      enddo
      enddo
      
      do i=ish+1,nx-ish
      do j=0,ish-1
      do m=1,ny
      x(i,m)=x(i,m)+sh(j)*(y(i+ish-j,m)+y(i-ish+j,m))
      enddo
      enddo
      do m=1,ny
      x(i,m)=(1.0+sh(ish))*y(i,m)+x(i,m)
      enddo
      enddo
      
      do i=nx-ish+1,nx-1
      do j=0,ish-1
      mm=(i+ish-j)-nx
      if(mm.GE.0)then
      do m=1,ny
      x(i,m)=x(i,m)+sh(j)*(2.0*y(nx,m)-y(nx-mm,m)+y(i-ish+j,m))
      end do
      else
      do m=1,ny
      x(i,m)=x(i,m)+sh(j)*(y(i+ish-j,m)+y(i-ish+j,m))
      end do
      endif
      enddo
      do m=1,ny
      x(i,m)=(1.0+sh(ish))*y(i,m)+x(i,m)
      end do
      enddo
      
      return
      end
      
