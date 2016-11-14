program bodc
implicit none
integer, parameter :: nx=1321
integer, parameter :: ny=1501
real bat(nx,ny),tmp(nx,ny)
real lonmin,lonmax,latmin,latmax,dlon,dlat
integer i,j

open(10,file='ammp_bathy.dat')
read(10,'(4f6.1,2i5)')lonmin,lonmax,latmin,latmax,i,j
print *,lonmin,lonmax,latmin,latmax,i,j
read(10,'(f6.1,10f7.1)')tmp(1:nx,1:ny)
close(10)
do j=1,ny
   bat(:,j)=tmp(:,ny-j+1)
enddo


dlon=(lonmax-lonmin)/float(nx-1); print '(a,2f8.5)','dlon=',dlon,dlon*60.0
dlat=(latmax-latmin)/float(ny-1); print '(a,2f8.5)','dlat=',dlat,dlat*60.0

open(10,file='tecbodc.dat')
   write(10,*)'TITLE = "BODC Bathymetry"'
   write(10,*)'VARIABLES = "i-index" "j-index" "Longitude" "Lattitude" "Bathymetry"'
   write(10,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',nx,', J=',ny,', K=1'

! Positions in lon-lat
   write(10,'(20I5)')((i,i=1,nx),j=1,ny)
   write(10,'(20I5)')((j,i=1,nx),j=1,ny)
   write(10,900)((lonmin+float(i-1)*dlon,i=1,nx),j=1,ny)
   write(10,900)((latmin+float(j-1)*dlat,i=1,nx),j=1,ny)
   write(10,900)((bat(i,j),i=1,nx),j=1,ny)
   900 format(10(1x,e12.5))
close(10)

end program bodc

