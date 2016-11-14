module m_nearestpoint
contains
subroutine nearestpoint(glon,glat,nx,ny,lon,lat,ipiv,jpiv,a1,a2,a3,a4,gmsk,ass)
   implicit none

   integer, intent(in)  :: nx,ny
   real,    intent(in)  :: glon(nx,ny),glat(nx,ny)
   real,    intent(in)     :: lon,lat
   integer, intent(inout)  :: ipiv,jpiv
   real,    intent(out)    :: a1,a2,a3,a4
   logical, intent(in)     :: gmsk(nx,ny)
   logical, intent(out)    :: ass

   real, allocatable  :: A(:,:)
   integer isize
   integer ia,ib,ja,jb
   integer iloc(2),i,j
   real, external :: spherdist

   iloc=-1

   do isize=1,20,2
      ia=max(1,ipiv-isize)
      ib=min(nx,ipiv+isize)
      ja=max(1,jpiv-isize)
      jb=min(ny,jpiv+isize)

      allocate( A(ia:ib,ja:jb) ) 
      
      do j=ja,jb
      do i=ia,ib
         if (gmsk(i,j)) then
            A(i,j)=spherdist(glon(i,j),glat(i,j),lon,lat)
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

!   print *,'OLD pivots:',ipiv,jpiv
   ipiv=ia+iloc(1)-1
   jpiv=ja+iloc(2)-1
!   print *,'NEW pivots:',ipiv,jpiv
   if (gmsk(ipiv,jpiv)) then
      a1=1.0
      ass=.true.
   else
      print *,'nearest_point (WARNING): Could not find (in,jn) for ipiv, jpiv=',ipiv,jpiv
      a1=0.0
      ass=.false.
   endif
   
   a2=0.0
   a3=0.0
   a4=0.0
end subroutine nearestpoint
end module m_nearestpoint
