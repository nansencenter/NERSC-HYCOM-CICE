module m_topofix
contains
subroutine topofix(depths,nx,ny)

   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   real, intent(inout) :: depths(nx,ny)
   integer firsti,firstj
   integer lasti,lastj
   logical :: ex
   integer :: ios
   real    :: depth
   character(len=80) :: c80

   inquire(file='grid.topofix',exist=ex)
   if (ex) then
      open(10,file='grid.topofix',status='old')
      print *,'Applying topofix corrections'
      ios=0
      do while (ios==0)

         ! Check if line is "commented"
         read(10,'(a)',iostat=ios) c80
         if (c80(1:1)=='#') then
            print *,'Skipping this line :'//trim(c80)
            cycle ! Skip to next execution of while loop
         else
            backspace(10) ! this shuld contain real data, skip back and read below
         end if

         read(10,*,iostat=ios) firsti,lasti, firstj,lastj,depth
         if (ios==0) then
            depths(firsti:lasti,firstj:lastj)=depth
            print '("depths(",i5,":",i5,",",i5,":",i5,")=",f7.1)',&
               firsti,lasti,firstj,lastj,depth
         end if
      end do
      close(10)
   else
      print *,'No topofix corrections'
   end if


end subroutine topofix
end module m_topofix
