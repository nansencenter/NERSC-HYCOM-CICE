module m_nearestpoint
contains
subroutine nearestpoint(glon,glat,nx,ny,lon,lat,ipiv,jpiv,a1,a2,a3,a4,gmsk,ass,lperiodic)
   implicit none

   integer, intent(in)  :: nx,ny
   real,    intent(in)  :: glon(nx,ny),glat(nx,ny)
   real,    intent(in)     :: lon,lat
   integer, intent(inout)  :: ipiv,jpiv
   real,    intent(out)    :: a1,a2,a3,a4
   logical, intent(in)     :: gmsk(nx,ny)
   logical, intent(out)    :: ass
   logical, intent(in )    :: lperiodic

   real, allocatable  :: A(:,:)
   integer isize
   integer ia,ib,ja,jb,ipiv0,jpiv0
   integer iloc(2),i,j,newi
   real, external :: spherdist

   ipiv0=ipiv
   jpiv0=jpiv
   iloc=-1
   do isize=1,20,2
      if (.not.lperiodic) then
         ia=max(1,ipiv-isize)
         ib=min(nx,ipiv+isize)
      else ! periodic in i
         ia=ipiv-isize
         ib=ipiv+isize
      end if
      ja=max(1,jpiv-isize)
      jb=min(ny,jpiv+isize)

      allocate( A(ia:ib,ja:jb) ) 
      
      do j=ja,jb
      do i=ia,ib
         if (lperiodic) then
            newi=mod(5*nx+i-1,nx)+1
         else
            newi=i
         end if
         if (gmsk(newi,j)) then
            A(i,j)=spherdist(glon(newi,j),glat(newi,j),lon,lat)
         else
            A(i,j)=1.0E20
         endif
      enddo
      enddo

      if (minval(A) < 1.0E20) then
         iloc=minloc(A)
         deallocate(A)
         exit
      endif
      deallocate(A)
   enddo

   !print *,'OLD pivots:',ipiv,jpiv
   if (lperiodic) then
      ipiv=mod(ia+iloc(1)-1+nx-1,nx)+1
   else
      ipiv=ia+iloc(1)-1
   end if
   jpiv=ja+iloc(2)-1
   !print '(a,6i5)','NEW pivots:',ipiv,jpiv,iloc,ia,ja
   if (all(iloc==-1)) then
      !print '(a,2i5)','(1)nearest_point (WARNING): Could not find (in,jn) for ipiv, jpiv=',ipiv0,jpiv0
      a1=0.0
      ass=.false.
   else if (gmsk(ipiv,jpiv)) then
      a1=1.0
      ass=.true.
   !else
   !   !print '(a,2i5)','(2)nearest_point (WARNING): Could not find (in,jn) for ipiv, jpiv=',ipiv0,jpiv0
   !   a1=0.0
   !   ass=.false.
   endif
   
   a2=0.0
   a3=0.0
   a4=0.0
end subroutine nearestpoint
end module m_nearestpoint
