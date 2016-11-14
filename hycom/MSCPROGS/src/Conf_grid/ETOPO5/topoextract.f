program topography
   implicit none

   integer, parameter :: nx=181     ! i-dim of model grid
   integer, parameter :: ny=240     ! j-dim of model grid

   integer, parameter :: nrx=4320   ! xdim of ETOPO5
   integer, parameter :: nry=2160   ! ydim of ETOPO5
   integer ds759(nrx,nry)           ! ETOPO5 file
   integer etopo5(nrx,nry)           ! ETOPO5 file

   real rad
   integer i,j,ia,ib,ja,jb
   real lonmin,lonmax,latmin,latmax
   real, allocatable, dimension(:,:) :: lon,lat,x,y,z

   rad=4.*atan(1.)/180.
   open(10, file='infile')
      read(10,*)lonmin
      read(10,*)lonmax
      read(10,*)latmin
      read(10,*)latmax
   close(10)
        
   ia=nint(lonmin*12.0+1.0)
   ib=nint(lonmax*12.0+1.0)
   jb=nint((90.0-latmin)*12.0+1.0)
   ja=nint((90.0-latmax)*12.0+1.0)
   if (ib==nrx+1) ib=nrx
   if (jb==nry+1) jb=nry

   print *,ia,ib,ja,jb

   allocate(lon(ia:ib,ja:jb))
   allocate(lat(ia:ib,ja:jb))
   allocate(x(ia:ib,ja:jb))
   allocate(y(ia:ib,ja:jb))
   allocate(z(ia:ib,ja:jb))

! Reading DS759.2 data set
   open(10,file='DS759.2.uf',form='unformatted')
      read(10)ds759
   close(10)
   write(*,*)'DS759.2 is read'

! Reading ETOPO5 data set
   open(10,file='ETOPO5.uf',form='unformatted')
      read(10)etopo5
   close(10)
   write(*,*)'ETOPO5 is read'

   do j=ja,jb
   do i=ia,ib
      lon(i,j)=float(i-1)/12.0
      lat(i,j)=90.0-float(j-2)/12.0
      x(i,j)=cos(lat(i,j)*rad)*cos(lon(i,j)*rad)
      y(i,j)=cos(lat(i,j)*rad)*sin(lon(i,j)*rad)
      z(i,j)=sin(lat(i,j)*rad)
      etopo5(i,j)= abs(min(etopo5(i,j),0))
      ds759(i,j)= abs(min(ds759(i,j),0))
   enddo
   enddo

   open(10,file='topo.dat')
      write(10,*)'TITLE = "TOPO"'
      write(10,*)'VARIABLES = "i"  "j"  "lon" "lat" "x" "y" "z" "ds759" "etopo5 "diff"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ib-ia+1,', J=',jb-ja+1,', K=1'

! Positions in Grid index
      write(10,'(10I8)')((i,i=ia,ib),j=ja,jb)
      write(10,'(10I8)')((j,i=ia,ib),j=ja,jb)

! Positions in Lon Lat
      write(10,900)((lon(i,j),i=ia,ib),j=ja,jb)
      write(10,900)((lat(i,j),i=ia,ib),j=ja,jb)
      write(10,900)((x(i,j),i=ia,ib),j=ja,jb)
      write(10,900)((y(i,j),i=ia,ib),j=ja,jb)
      write(10,900)((z(i,j),i=ia,ib),j=ja,jb)
      write(10,'(10I8)')((ds759(i,j),i=ia,ib),j=ja,jb)
      write(10,'(10I8)')((etopo5(i,j),i=ia,ib),j=ja,jb)
      write(10,'(10I8)')((ds759(i,j)-etopo5(i,j),i=ia,ib),j=ja,jb)
   close(10)

900  format(10(1x,e12.5))

end program
