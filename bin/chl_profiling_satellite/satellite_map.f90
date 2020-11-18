! to compile execute the following in the command line:
! f2py -c -m _the_mapping_loop_with_ice satellite_map.f90 
subroutine main(bathy,scpx,scpy,plat,plon,icemask,&
           slat,slon,chlsat,satint,II,JJ,IIs,JJs)
   implicit none
   integer, intent(in) :: II,JJ,IIs,JJs
   real, intent(out),dimension(JJ,II) :: satint
   real, intent(in),dimension(JJ,II) :: plat,plon,scpx,scpy,bathy
   real, intent(in),dimension(JJs) :: slat
   real, intent(in),dimension(IIs) :: slon
   real, intent(in),dimension(JJs,IIs) :: chlsat
   integer, intent(in),dimension(JJ,II) :: icemask
   integer :: I,J,indlat,indlon,search,counted,IM,JM
   real :: chlmean,arc
!f2py intent(in) :: II,JJ,KK,IIs,JJs
!f2py intent(in) :: plat,plon,slat,slon,chlsat,bathy,icemask
!f2py intent(out) :: satint

satint    = 1E36
DO IM=1,II
  DO JM=1,JJ
    indlat = minloc(abs(slat - plat(JM,IM)),dim=1)
    indlon = minloc(abs(slon - plon(JM,IM)),dim=1)
    ! A small satellite search region around the JM,IM point is defined
    ! grid size/8000 corresponds to (grid size/2 - half distance in the grid)/(4km satellite res.)
    ! multiplication by 2 at the end extends the search area to capture enough
    ! points 
    search = nint(max(scpx(JM,IM),scpy(JM,IM))/8000.) * 2
    counted = 0
    chlmean = 0
    IF (bathy(JM,IM)>=100. .and. bathy(JM,IM)<2E5 .and. icemask(JM,IM)==0) THEN 
    ! process model points deeper than 100m and no ice
      DO I=max(indlon-search,1),min(indlon+search,IIs)
         DO J=max(indlat-search,1),min(indlat+search,JJs)

          call distance(slat(J),plat(JM,IM),slon(I),plon(JM,IM),arc)
          IF (arc <= max(scpx(JM,IM),scpy(JM,IM))/2000.) THEN
            IF (chlsat(J,I)<1E4 .and. chlsat(J,I)>=0.) THEN
                  chlmean = chlmean + chlsat(J,I)
                  counted = counted + 1
            ENDIF
          ENDIF
         ENDDO
      ENDDO
      IF (counted ==0) THEN
         chlmean = 1E36
      ELSE
         chlmean = chlmean / counted
      ENDIF
      IF (chlmean==0.) chlmean = 1E36

      satint(JM,IM) = chlmean

    ELSE ! bathymetry check
      satint(JM,IM)    = 1E36 ! this high value will be masked out by python
                              ! so no profiling at these points
    ENDIF ! bathymetry check
  ENDDO
ENDDO
end subroutine main




subroutine distance(lat1,lat2,long1,long2,arc)
   implicit none
   real, parameter :: deg2rad = 3.1415927 / 180.
   real :: phi1,phi2,theta1,theta2   
   real, intent(in)  :: lat1,lat2,long1,long2
   real, intent(out) :: arc 

   phi1 = (90.0 - lat1)*deg2rad
   phi2 = (90.0 - lat2)*deg2rad

   theta1 = long1*deg2rad
   theta2 = long2*deg2rad

   arc =  acos( sin(phi1) * sin(phi2) * cos(theta1 - theta2) + &
          cos(phi1) * cos(phi2) )
   arc = arc * 6378.137 ! convert to km

end subroutine distance 
