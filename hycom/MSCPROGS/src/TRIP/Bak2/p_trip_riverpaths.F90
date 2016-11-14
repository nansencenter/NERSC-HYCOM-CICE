! --- -------------------------------------------------------------------
! --- Diagnostic routine trip_riverpaths
! --- -------------------------------------------------------------------
! --- Program to "connect the dots" in the trip05 river path database.
! --- For each land point, the trip database is used to find where water
! --- originating here ends up in the ocean. This routine is mainly used
! --- for diagnostics - although the catchment table can be used by the 
! --- netcdf output in trip_flow
! ---
! --- Output from this routine is:
! --- land point location, final destination, basin flag
! --- -------------------------------------------------------------------

program trip_riverpaths
   use mod_trip
   implicit none
   integer, parameter :: maxniter=2*(nx+ny) ! Should do
   real, parameter :: rearth=6372.795477598 ! Quadratic mean radius (km)
   real, parameter :: radian=57.2957795

   real   , dimension(nx,ny) :: distance
   integer, dimension(nx,ny) :: destinationx, destinationy,lmask
   real :: tmplon,tmplat
   integer :: tmpdir, ios
   integer :: i,j,i2,j2,i3,j3, niter
   logical :: same, ex
   integer :: basins(nx,ny), dirx(nx,ny), diry(nx,ny)
   integer :: basinx(nx*ny), basiny(nx*ny), sumbasin(nx*ny)
   real    :: sumcatch(nx*ny),tmp
   integer :: nbasin, ibasin, ibasin2
   character(len=200) :: cenv,path0


   ! set trip
   call init_trip()

   ! TRIP grid
   ! TODO: This should be done in trip_init - but for now there
   !       is som inconsistency between "flips" in riverweights/riverflow
   do j=1,ny
      lat(j) = 90-dx/2 - (j-1)*dy 
   end do
   do i=1,nx
      lon(i) = (i-1)*dx + dx/2
   end do


   ! Start connecting the dots - 
   print *,'Connecting TRIP grid cells with Sea cells'
   destinationx=-1
   destinationy=-1
   dirx=0
   diry=0
   distance=0.
   do j=1,ny
   do i=1,nx
   if (direction(i,j)/=9.and.direction(i,j)/=0) then ! land point with throughflow

      i2=i
      j2=j
      i3=-1
      j3=-1
      niter=0

      do while (direction(i2,j2)/=9.and.direction(i2,j2)/=0.and.niter<=maxniter)

         i3=i2
         j3=j2
         if     (direction(i3,j3)==1.or.direction(i3,j3)==2.or.  direction(i3,j3)==8) then
            j2=max(j2-1,1)
            diry(i3,j3)=1
         elseif (direction(i3,j3)==4.or.direction(i3,j3)==5.or.  direction(i3,j3)==6) then
            j2=min(j2+1,ny)
            diry(i3,j3)=-1
         end if


         if     (direction(i3,j3)==2.or.direction(i3,j3)==3.or.  direction(i3,j3)==4) then
            i2=mod(i2,nx)+1
            dirx(i3,j3)=1
         elseif (direction(i3,j3)==6.or.direction(i3,j3)==7.or.  direction(i3,j3)==8) then
            i2=mod(nx+i2-2,nx)+1
            dirx(i3,j3)=-1
         end if

         distance(i,j)    = distance(i,j)  +                     &
                            sqrt(real( max(-1,min(i2-i3,1))**2 + &
                                       max(-1,min(j2-j3,1))**2 ))



         ! Safety check
         !if (niter==maxniter) then
         !   print *,'too many iterations in traversion of grid '
         !   print *,i,j,i2,j2,niter,maxniter
         !!   !stop '(trip_riverpaths)'
         !end if

         niter=niter+1
         !print *,i2,j2
      end do

      if (niter<maxniter) then
         destinationx(i,j)=i2
         destinationy(i,j)=j2
      else
         destinationx(i,j)=-1
         destinationy(i,j)=-1
         distance(i,j)    = -100
      end if


   end if
   end do
   end do


   ! Find unique basins
   nbasin=0
   basins=-1
   sumbasin=0
   sumcatch=0
   do j=1,ny
   do i=1,nx
      if (destinationx(i,j)/=-1) then

         ibasin2=-1
         do ibasin=1,nbasin
            if (destinationx(i,j)==basinx(ibasin) .and.  &
                destinationy(i,j)==basiny(ibasin)) then
               ibasin2=ibasin
            end if
         end do


         ! New basin
         if (ibasin2==-1) then
            nbasin=nbasin+1
            basinx(nbasin)=destinationx(i,j)
            basiny(nbasin)=destinationy(i,j)
            sumbasin(nbasin)=1
            basins(i,j)=nbasin
         else
            basins(i,j)=ibasin2

            ! sum of grid points ending up here
            sumbasin(ibasin2)=sumbasin(ibasin2)+1

            ! Catchment area
            sumcatch(ibasin2)=sumcatch(ibasin2) + &
               abs(cos(lat(j)/radian))*cos(dx/radian)*cos(dy/radian)*rearth**2
         end if
      end if
   end do
   end do

   print '(a,i6)','nbasins',nbasin !,count(destinationx/=-1)


   print '(a)','trips.tec        - table of water origin / water exit grid cells'
   open(10,file='trips.tec',status='replace')
   do j=1,ny
   do i=1,nx
      write(10,'(2i5,2e14.4,2i5,e14.4,4i5)') i,j,lon(i),lat(j), &
         destinationx(i,j),destinationy(i,j),distance(i,j),direction(i,j), &
         basins(i,j), dirx(i,j), diry(i,j)
   end do
   end do
   close(10)


      
   print '(a)','trips_rivers.tec - table of river outlet  catchment area'
   open(10,file='trips_rivers.tec',status='replace')
   do i=1,nbasin
      write(10,'(2i5,2e14.4,2i5,e14.4)') basinx(i),basiny(i), &
         lon(basinx(i)),lat(basiny(i)),i,sumbasin(i), sumcatch(i)
   end do
   close(10)

   ! Sort basins on catchment area (slow sort yes)
   do i =1,nbasin
   do i2=i,nbasin
   if (sumcatch(i2)>sumcatch(i)) then
      tmp=basinx  (i) ; basinx  (i)=basinx  (i2); basinx  (i2)=tmp
      tmp=basiny  (i) ; basiny  (i)=basiny  (i2); basiny  (i2)=tmp
      tmp=sumcatch(i) ; sumcatch(i)=sumcatch(i2); sumcatch(i2)=tmp
   end if
   end do
   end do




   !Top rivers by catchment
   print '(a)','catchment.asc    - top rivers by catchment'
   open(10,file='catchment.asc',status='replace')
   do i=1,nbasin
      write(10,'(2i5,2f14.5,e14.4)') basinx(i),basiny(i), &
         lon(basinx(i)),lat(basiny(i)), sumcatch(i)
   end do
   close(10)






   end program
