
subroutine lon_lat_extrap(modlat,modlon,plat,plon,nx,ny)
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   real, intent(in) :: modlat(nx,ny)
   real, intent(in) :: modlon(nx,ny)
   real, intent(out) :: plat(0:nx+1,0:ny+1)
   real, intent(out) :: plon(0:nx+1,0:ny+1)

   plat(1:nx,1:ny)=modlat(1:nx,1:ny)
   plon(1:nx,1:ny)=modlon(1:nx,1:ny)

   plat(0,:)=4.0*plat(1,:)-6.0*plat(2,:)+4.0*plat(3,:)-plat(4,:)
   plat(nx+1,:)=4.0*plat(nx,:)-6.0*plat(nx-1,:)+4.0*plat(nx-2,:)-plat(nx-3,:)
   plat(:,0)=4.0*plat(:,1)-6.0*plat(:,2)+4.0*plat(:,3)-plat(:,4)
   plat(:,ny+1)=4.0*plat(:,ny)-6.0*plat(:,ny-1)+4.0*plat(:,ny-2)-plat(:,ny-3)

   plon(0,:)=4.0*plon(1,:)-6.0*plon(2,:)+4.0*plon(3,:)-plon(4,:)
   plon(nx+1,:)=4.0*plon(nx,:)-6.0*plon(nx-1,:)+4.0*plon(nx-2,:)-plon(nx-3,:)
   plon(:,0)=4.0*plon(:,1)-6.0*plon(:,2)+4.0*plon(:,3)-plon(:,4)
   plon(:,ny+1)=4.0*plon(:,ny)-6.0*plon(:,ny-1)+4.0*plon(:,ny-2)-plon(:,ny-3)
end subroutine lon_lat_extrap
