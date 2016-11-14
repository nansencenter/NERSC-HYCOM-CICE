module m_bigrid
contains
subroutine bigrid(depths,nx,ny)
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   real, intent(inout) :: depths(nx,ny)
   logical periodic
   integer nfill,nzero,jsec,jfrst,jlast,jjnew,ii,ii1,jj,jj1
   integer l,j,i,ja,jb
   ii1=nx-1
   jj1=ny-1
   jj=ny
   ii=nx

   print *,'Bigrid iterations:'
   do l=1,100
      nfill=0
      !write(*,*)'bigrid iteration: ',l

      ! KAL(1) do j=1,jj1 
      do j=1,jj !KAL(1)
         ja=mod(j-2+jj,jj)+1
         jb=mod(j     ,jj)+1
         ! KAL(1) do i=1,ii1
         do i=1,ii !KAL(1)
            nzero=0
            if (depths(i,j).gt.0.) then
               if (i.eq.1  .or.depths(max(1,i-1),j).le.0.) nzero=nzero+1
               if (i.eq.ii1.or.depths(i+1,j).le.0.) nzero=nzero+1
               ! KAL(1) if (j.eq.1  .or.depths(i,ja).le.0.) nzero=nzero+1 
               ! KAL(1) if (j.eq.jj1.or.depths(i,jb).le.0.) nzero=nzero+1
               if (depths(i,ja).le.0.) nzero=nzero+1 !KAL(1)
               if (depths(i,jb).le.0.) nzero=nzero+1 !KAL(1)
               if (nzero >= 3) then
                  !write (*,'(a,i4,a,i4,a)')' depths(',i,',',j,') set to zero'
                  depths(i,j)=0.
                  nfill=nfill+1
               endif
            endif
         enddo
      enddo
      print '("Iteration",i4," filled ",i7," points")', l,nfill
      if (nfill == 0) exit
   enddo
   !KAL(1) - Changes to allow for periodic grid (in j-direction) -- 25.10.2005
end subroutine bigrid
end module m_bigrid
