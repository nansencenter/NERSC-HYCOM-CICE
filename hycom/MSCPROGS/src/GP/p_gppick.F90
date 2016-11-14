program gppick
   implicit none
   integer i,mini
   character(len=80) gpdir
   real lon,lat,gplon,gplat
   integer nrgp
   character(len=3) cstation
   character(len=1) sgn
   real dist
   real mindist
   character(len=80), allocatable :: gpfile(:)
   logical ex
   real, external :: spherdist

   inquire(exist=ex,file='inp.gp')
   if (ex) then
      open(10,file='inp.gp')
      read(10,'(a)')gpdir
      read(10,*)lon
      read(10,*)lat
      read(10,*)cstation
      close(10)
   else
      write(*,'(a)',advance='no')'Choose GP directory: '
      read(*,*)gpdir
      write(*,'(a)',advance='no')'Choose Longitude : '
      read(*,*)lon
      write(*,'(a)',advance='no')'Choose Latitude  : '
      read(*,*)lat
      write(*,'(a)',advance='no')'"cstation" (station id): '
      read(*,*) cstation
   end if

   call system('ls '//trim(gpdir)//' | grep "gp_" > tmp_gp')
   call system('cat tmp_gp | wc -l  > tmp_gpnr')

   open(10,file='tmp_gpnr')
      read(10,*)nrgp
   close(10)


   print *,'number of files are: ',nrgp
   allocate(gpfile(nrgp))

   open(10,file='tmp_gp')
   open(11,file='tmp_gp.out')

   mindist=100000000.0
   do i=1,nrgp
      read(10,'(a80)',err=100,end=200)gpfile(i)
      !print *,gpfile(i)
      read(gpfile(i)(12:12),'(a1)',err=100)sgn
      read(gpfile(i)(13:19),'(f7.3)',err=100)gplon
      if (sgn=='-') gplon=-gplon
      read(gpfile(i)(21:21),'(a1)',err=100)sgn
      read(gpfile(i)(22:28),'(f7.3)',err=100)gplat
      if (sgn=='-') gplat=-gplat
      !print *,i,gplon,gplat
      dist=spherdist(lon,lat,gplon,gplat)
      write(11,'(f13.2,4f9.3,tr2,a)')dist*1e-3,lon,lat,gplon,gplat,trim(gpfile(i))
      if (dist < mindist) then
         mindist=dist
         mini=i
      endif
      100 continue
   enddo

   200 close(10)
   close(11)

   ! Clean up
   call system('rm tmp_gp')
   call system('rm tmp_gpnr')

   call system('cat tmp_gp.out | sort > tmp_gp.sort')
   print *,cstation,' ',trim(gpfile(mini)),' ',mindist*1e-3
   print *,'pos info in filename: ',gpfile(mini)(12:28)

   inquire(file='infofile.txt',exist=ex)
   if (.not.ex) then
      open(10,file='infofile.txt',status='new')
         write(10,'(a3,5(a10))')' NO','   GPP-LON','   GPP-LAT','   CPP-LON','   CPP-LAT','  DIST(km)'
      close(10)
   endif
   open(10,file='infofile.txt',position='append')
      write(10,'(a3,2(tr2,a8),2(tr2,f8.3),f12.2)')cstation,gpfile(mini)(12:19), &
         gpfile(mini)(21:28),lon,lat,mindist*1e-3
   close(10)
end program






