program gppick2
! Alternate version of gppick - based on file contents rather than 
! name
   use mod_read_gp
   implicit none
   integer i,mini
   character(len=80) gpdir
   real lon,lat,gplon,gplat
   integer nrgp, ios, irec
   character(len=3) cstation, rungen
   real dist
   real mindist, minlon, minlat
   character(len=80), allocatable :: gpfile(:)
   character(len=80) :: tmparg
   logical ex, okrec
#if defined (IARGC)
   integer*4, external :: iargc
#endif
   real, external :: spherdist

   if (iargc()==5) then
      call getarg(1,rungen);
      call getarg(2,tmparg); read(tmparg,*) lon
      call getarg(3,tmparg); read(tmparg,*) lat
      call getarg(4,cstation);
      call getarg(5,gpdir)
   else
      print *
      print *,'gppick2 searches through files in <directory> to get the '
      print *,'GP station closest to the specified lon/lat position     '
      print *,'also needed is a station ID giving the name of the station'
      print *,'This station name is up to you!'
      print *
      print *,'Usage: gppick2 <rungen> <lon> <lat> <station id> <directory>'
      print *
      print *,'The file <rungen>gprecord.asc must be present in <directory>'
      call exit(1)
   end if


   ! Get info in record
   call read_gpheader(rungen,trim(gpdir))

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

      if (mod(i,100)==0) print'(a,i5,a,i6)','processed ',i,' of ',nrgp

      ! File name to be checked
      read(10,'(a80)',err=100,end=200)gpfile(i)(:)
      irec=1
      ios=0
      okrec=.false.
      do while (ios==0)
         call read_gprecord(trim(gpdir)//'/'//trim(gpfile(i)),irec,ios)
         if (ios==0.and.itime>0) then
            gplon  =readvar('lon')
            gplat  =readvar('lat')
            okrec=.true.
            exit
         end if
         irec=irec+24
      end do

      ! If the read above failed some how
      if (.not.okrec) cycle

      dist=spherdist(lon,lat,gplon,gplat)
      write(11,'(f13.2,4f9.3,tr2,a)')dist*1e-3,lon,lat,gplon,gplat,trim(gpfile(i))
      if (dist < mindist) then
         mindist=dist
         mini=i
         minlon=gplon
         minlat=gplat
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
   print *,'pos info in data file:',minlon,minlat

   inquire(file='infofile.txt',exist=ex)
   if (.not.ex) then
      open(10,file='infofile.txt',status='new')
      write(10,'(a3,6(a10))')' NO','   GPP-LON','   GPP-LAT','   CPP-LON','   CPP-LAT','  DIST(km)','      Name'
      close(10)
   endif
   open(10,file='infofile.txt',position='append')
      write(10,'(a3,2(tr2,a8),2(tr2,f8.3),f10.2,tr2,a)')cstation,gpfile(mini)(12:19), &
         gpfile(mini)(21:28),lon,lat,mindist*1e-3,trim(gpfile(mini))
   close(10)
end program






