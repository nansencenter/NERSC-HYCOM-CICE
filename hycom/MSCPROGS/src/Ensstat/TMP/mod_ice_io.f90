module mod_ice_io



contains


! Slim version of the routine in used in HYCOM
   subroutine read_restart_ice(fname,fld,cfld,nx,ny,imem)
   implicit none
   integer,                intent(in)  :: nx,ny,imem
   character(len=*),       intent(in)  :: cfld,fname
   real, dimension(nx,ny), intent(out) :: fld

   real*8, dimension(nx,ny) :: &
     ficem,hicem,hsnwm,ticem,tsrfm

   logical :: ex
   integer j,iostat
      

   inquire(file=fname,exist=ex)

   if (ex) then
      inquire(iolength=j) &
         ficem,hicem,hsnwm,ticem,tsrfm  
      open(10,file=fname,status='old',form='unformatted', &
              access='direct',recl=j)
      read(10,rec=imem,err=200,iostat=iostat)ficem,hicem, &
                 hsnwm,ticem,tsrfm 
  200 continue
      close(10)
      if (iostat/=0) then
         write(6,*) 'An error occured while reading '//trim(fname)
         write(6,*) 'IOSTAT is :',iostat
      else

         if (trim(cfld)=='hicem') then
            fld=hicem
         elseif (trim(cfld)=='ficem') then
            fld=ficem
         elseif (trim(cfld)=='hsnwm') then
            fld=hsnwm
         else
            print *,'Unknown field '//trim(cfld)
            stop '(read_restart_ice)'
         end if
      end if

   else
      write(6,*) 'File not found ...'//trim(fname)
      stop '(read_restart_ice)'
   endif
   end subroutine read_restart_ice



! Reads ice drift fields, dumped by daily average module
   subroutine read_icedrift(fname,fld,cfld,nx,ny,imem,undef)
   implicit none
   integer,                intent(in)  :: nx,ny,imem
   character(len=*),       intent(in)  :: cfld,fname
   real, dimension(nx,ny), intent(out) :: fld
   real                  , intent(in)  :: undef

   real*4, dimension(nx,ny) :: idrftu,idrftv

   logical :: ex
   integer j,iostat
      

   inquire(file=fname,exist=ex)

   if (ex) then
      inquire(iolength=j) idrftu, idrftv
         
      open(10,file=fname,status='old',form='unformatted', &
              access='direct',recl=j)
      read(10,rec=imem,err=200,iostat=iostat) idrftu, idrftv
  200 continue
      print *,'read_icedrift iostat',iostat
      print *,'read_icedrift iostat',iostat
      print *,'read_icedrift iostat',iostat
      print *,'read_icedrift iostat',iostat
      print *,'read_icedrift iostat',iostat
      close(10)
      if (iostat/=0) then
         write(6,*) 'An error occured while reading '//trim(fname)
         write(6,*) 'IOSTAT is :',iostat
      else

         if (trim(cfld)=='drfticex') then
            fld=idrftu
         elseif (trim(cfld)=='drfticey') then
            fld=idrftv
         else
            print *,'Unknown field '//trim(cfld)
            stop '(read_icedrift)'
         end if

         write(6,*) 'Got field '//trim(cfld)
         where(abs(fld-undef)/abs(undef)<1e-6) fld=0.
      end if

   else
      write(6,*) 'File not found ...'//trim(fname)
      stop '(read_icedrift)'
   endif
   end subroutine read_icedrift






end module mod_ice_io
