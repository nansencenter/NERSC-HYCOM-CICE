module m_checktopo
contains
subroutine checktopo(depths,nx,ny)

   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   real, intent(in)    :: depths(nx,ny)
   logical lmsk(nx,ny),donep(nx,ny),donetmp(nx,ny)
   integer i,j,ia,ja,iia,jja,l,isign
   integer ifalse, ifalseold

   lmsk=.false.
   donep=.false.
   where (depths > 0.0) lmsk=.true.
   j=ny/2
   do i=1,nx
      if (depths(i,j) /= 0.0) then
         lmsk(i,j)=.true.
         exit
      endif
   enddo

      ifalseold=0
   do l=1,1000
      if (mod(l,2) == 0) isign=1
      if (mod(l,2) == 1) isign=-1

      ifalse=0
      donetmp=.false.
      do j=1,ny
      do i=1,nx
         if (lmsk(i,j).and.(.not.donep(i,j))) then
            ifalse=ifalse+1
            do jja=-1,1
               ja=j+jja*isign
               if (ja < 1) ja=ny
               if (ja > ny) ja=1

               do iia=-1,1
                  ia=i+iia*isign
                  if (ia < 1) ia=1
                  if (ia > nx) ia=nx

                  if (lmsk(ia,ja).and.donep(ia,ja)) then
                     donetmp(i,j)=.true.
                     ifalse=ifalse-1
                     exit
                  endif
               enddo
               if (donetmp(i,j)) exit
            enddo
         endif
      enddo
      enddo
      where (donetmp) donep=.true.

      if (ifalse == 0) exit

      if (ifalse == ifalseold) exit
      ifalseold=ifalse

   enddo
   if (ifalse /= 0) then
      do j=1,ny
      do i=1,nx
         if (lmsk(i,j).and.(.not.donep(i,j))) then
            print '(a,2I6)',' Ifalse problem in grip point (i,j): ',i,j
         endif
      enddo
      enddo
   endif



end subroutine checktopo
end module m_checktopo
