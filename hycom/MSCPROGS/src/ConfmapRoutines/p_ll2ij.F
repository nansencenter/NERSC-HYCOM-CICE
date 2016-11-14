      program ll2ij
      use mod_confmap
      implicit none
      integer*4, external :: iargc
      real    :: lon, lat
      real    :: lon_n, lat_n
      integer :: ipiv, jpiv
      character(len=20) :: tmparg
      integer :: nx,ny

      if (iargc()/=2) then
         print *,'Routine converts from longitude/latitude to'
         print *,'grid indexes. Grid index printed to stdout '
         print *,'if on model grid, printed to stderr if not '
         print *,'   Usage : confmap_ll2ij longitude latitude'
         call exit(1)
      end if
      call getarg(1,tmparg) ; read(tmparg,*) lon
      call getarg(2,tmparg) ; read(tmparg,*) lat
c
c --- Initialize mapping
      call initconfmap(nx,ny,.true.)
c
c --- Go from lon lat to conformal mapping lon lat
      call oldtonew(lat,lon,lat_n,lon_n)
c
c --- Go from conformal mapping lon lat to grid pivot points
      call pivotp(lon_n,lat_n,ipiv,jpiv,0.)
c 
      if (ipiv<nx .and. ipiv >=1 .and. jpiv<ny .and.jpiv>=1) then
c ---    Tell STDOUT
         write(6,'(2i6)') ipiv,jpiv
      else
c ---    Tell STDERR
         write(0,'(a,2i6)') 'Outside grid : ',ipiv,jpiv
      end if
c
      end program





