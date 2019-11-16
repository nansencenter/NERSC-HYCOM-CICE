module m_grid_consistency 
contains 

   subroutine grid_consistency(depths,nxl,nyl,inest)
   implicit none

   integer , intent(in)    :: nxl,nyl,inest
   real    , intent(inout) :: depths(nxl,nyl)

   integer :: i2p1,j2p1,i2m1,j2m1,i2,j2,i,j
   integer :: ibnd, ibnd2, bdim, idir, bndcnt, iter, niter, radius
   integer :: maxradius, maxbnd
   integer, allocatable :: bndx(:),bndy(:)
   logical, allocatable :: connected(:,:), bndconnected(:)
   logical, dimension(4):: neighbourhood
   integer :: countiso,count3n






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check for grid consistency. Remove 3-land neighbour sea points, 
! and check for small isolated basins
   !integer, allocatable :: bndx(:),bndy(:)
   !logical, allocatable :: connected(:,:), bndconnected(:)
   maxradius= inest/2+1
   maxbnd   = 2*maxradius+1
   allocate(bndx(4*maxbnd))
   allocate(bndy(4*maxbnd))
   allocate(bndconnected(4*maxbnd))
   allocate(connected(nxl,nyl))
   
   ! Grid consistency iterations ...
   print *,'Checking for 3-land neighbour sea points and isolated basins'
   !print *,'NBNBNB test points inserted' ; !depths(266:270,13:15)=10. !; depths(245:266,14)=10.
   niter=5
   iter=0
   do while( iter<=niter .and. (count3n/=0 .or. countiso/=0))

      print *,'Iteration # ',iter+1

      ! Check for 3-neighbour points, fill them
      count3n=0
      do j=1,nyl
      do i=1,nxl
      if (depths(i,j)> .5) then
            
         neighbourhood=(/ depths(i  ,min(j+1,nyl))<.5, depths(i  ,max(j-1,1))<.5, &
                          depths(min(i+1,nxl),j  )<.5, depths(max(i-1,1),j  )<.5 /)
         if (count(neighbourhood)>=3) then
            print '(a,2i5,a)','   point ',i,j,' has +3 land neighbours'
            depths(i,j)=0.
            count3n=count3n+1
         end if
      end if
      end do
      end do




      ! Check for isolated areas within radius +- (inest/2 +1)
      countiso=0
      do j=2,nyl-1
      do i=2,nxl-1
      if (depths(i,j)>.5) then

         connected(i,j)=.true.
         do radius=1,inest/2+1

            ! Map boundary points to a vector
            bndcnt=0
            bdim  =2*radius+1
            do ibnd=-radius,radius
               bndcnt=bndcnt+1

               ! bottom
               bndx(bndcnt        ) = min(max(1,i+ibnd  ),nxl)
               bndy(bndcnt        ) = min(max(1,j-radius),nyl)

               ! right
               bndx(bndcnt +  bdim) = min(max(1,i+radius),nxl)
               bndy(bndcnt +  bdim) = min(max(1,j-ibnd  ),nyl)

               ! top    
               bndx(bndcnt +2*bdim) = min(max(1,i-ibnd  ),nxl)
               bndy(bndcnt +2*bdim) = min(max(1,j+radius),nyl)

               ! right
               bndx(bndcnt +3*bdim) = min(max(1,i-radius),nxl)
               bndy(bndcnt +3*bdim) = min(max(1,j+ibnd  ),nyl)
            end do

            ! Traverse boundary backwards and forwards along boundary 
            ! (idir) to get points connected to center
            do idir=0,1
            do ibnd=1,4*bdim
               
               if (idir==1) then
                  ibnd2=4*bdim - (ibnd-1)
               else
                  ibnd2=ibnd
               end if
               i2=bndx(ibnd2)
               j2=bndy(ibnd2)
               i2p1=min(i2+1,nxl)
               j2p1=min(j2+1,nyl)
               i2m1=max(1,i2-1)
               j2m1=max(1,j2-1)

               ! Set connected status for this cell
               if (depths(i2,j2)>.5) then
                  connected(i2,j2) =           &
                     connected(i2m1,j2  ) .or. &
                     connected(i2p1,j2  ) .or. & 
                     connected(i2  ,j2p1) .or. & 
                     connected(i2  ,j2m1) 
               end if

               ! Connectedness for boundary
               bndconnected(ibnd)=connected(i2,j2)
            end do
            end do
            !print *,'Connected points : ',count(bndconnected(1:4*bdim))
         end do
         if (count(bndconnected(1:4*bdim))==0) then
            print '(a,2i5,a)','   point ',i,j,' is isolated '
            depths(i,j)=0.
            countiso=countiso+1
         end if

         ! Set to unconnected for next run
         do i2=max(1,i-radius),min(nxl,i+radius)
         do j2=max(1,j-radius),min(nyl,j+radius)
            connected(i2,j2)=.false.
         end do
         end do
         
      end if
      end do
      end do
      iter=iter+1
   end do ! niter


   end subroutine grid_consistency
end module m_grid_consistency
