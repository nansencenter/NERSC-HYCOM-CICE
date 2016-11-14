module m_get_global_grid
contains
subroutine get_global_grid(lon,lat,depths)
   use mod_xc
   use mod_za ! zaiost must be called before running the sub
   implicit none

   character(len=*), parameter :: fdepthsa = 'regional.depth.a'
   character(len=*), parameter :: fdepthsb = 'regional.depth.b'
   character(len=*), parameter :: fgrida  = 'regional.grid.a'
   character(len=*), parameter :: fgridb  = 'regional.grid.b'

   real, intent(out),dimension(idm,jdm) :: lon,lat,depths

   character(len=7)  :: tag7
   character(len=80) :: c80
   logical :: exa, exb, exold
   integer :: t_idm, t_jdm, nop, i
   real*8, dimension(idm,jdm) :: io1,io2
   integer,   dimension(idm,jdm) :: mask
   real :: amax, bmax, amin, bmin

   mask=1

   inquire(exist=exa,file=fdepthsa)
   inquire(exist=exb,file=fdepthsb)
   write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
   inquire(file='gdepths'//tag7//'.uf',exist=exold)

   nop = 9
      
   ! First option, check for regional.depths file
   if (exa.and.exb) then
      open(10,file=fdepthsb)
      do i=1,5
         read(10,'(a80)') c80
         !print *,i,trim(c80)
      end do
      ! Read max/min depths
      read(10,'(a80)') c80
      read(c80(20:len_trim(c80)),*) bmin,bmax
      !print *,'2',bmin,bmax
      close(10)


      call zaiopf(fdepthsa,'old', nop)
      call zaiord(depths,mask,.false.,amin,amax,nop)
      if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
              abs(amax-bmax).gt.abs(bmax)*1.e-4     ) then
       write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
         'error - .a and .b files not consistent:', &
         '.a,.b min = ',amin,bmin,amin-bmin, &
         '.a,.b max = ',amax,bmax,amax-bmax
       stop '(get_global_grid)'
     endif
      !print *,'2',amin,amax
     call zaiocl(nop)
   else if (exold) then
      open (unit=10,file='gdepths'//tag7//'.uf',status='old', &
            form='unformatted')
      read(10)io1
      close(10)
      depths=io1
   else
      print *,'No depth files found'
      stop '(get_global_grid)'
   end if

   where (depths>1e28) depths=0.



    
   inquire(exist=exa,file=fgrida)
   inquire(exist=exb,file=fgridb)
   inquire(file='gnewpos.uf',exist=exold)

   if (exa .and. exb) then
      open(10,file=fgridb)
      call zaiopf(fgrida,'old', nop)
      read(10,*)
      read(10,*) t_idm
      read(10,*) t_jdm

      ! Longitude
      read(10,'(a80)') c80
      read(c80(20:len_trim(c80)),*) bmin,bmax

      call zaiord(lon,mask,.false.,amin,amax,nop)
      if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
              abs(amax-bmax).gt.abs(bmax)*1.e-4     ) then
       write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
         'error - .a and .b files not consistent:', &
         '.a,.b min = ',amin,bmin,amin-bmin, &
         '.a,.b max = ',amax,bmax,amax-bmax
       stop '(get_global_grid)'
      endif

      ! Latitude
      read(10,'(a80)') c80
      read(c80(20:len_trim(c80)),*) bmin,bmax
      close(10)

      call zaiord(lat,mask,.false.,amin,amax,nop)
      if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
              abs(amax-bmax).gt.abs(bmax)*1.e-4     ) then
       write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
         'error - .a and .b files not consistent:', &
         '.a,.b min = ',amin,bmin,amin-bmin, &
         '.a,.b max = ',amax,bmax,amax-bmax
       stop '(get_global_grid)'
      endif
      call zaiocl(nop)
   else if (exold) then

      open(10,file='gnewpos.uf',form='unformatted',status='old')
      read(10)io1,io2
      close(10)
      lat=io1
      lon=io2
      where (lon<0.) lon=lon+360.
   else
      print *,'No grid files found'
      stop '(get_global_grid)'
   end if
end subroutine get_global_grid
end module m_get_global_grid
