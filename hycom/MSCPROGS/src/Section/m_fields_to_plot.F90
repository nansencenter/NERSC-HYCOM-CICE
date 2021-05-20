module m_fields_to_plot


   type fields
      character(len=8) fextract
      integer          laya
      integer          layb
      logical          option
      logical          vecflag ! New - vector flag
      character(len=8) vecpost
   end type fields


contains
subroutine fields_to_plot(fld,nfld,hfile)
use mod_hycomfile_io
   implicit none
   integer, parameter  :: iversion=4
   type(fields), intent(out) :: fld(1000)
   integer, intent(out) :: nfld
   type(hycomfile), intent(in) :: hfile

   integer i,ivers,lay1,lay2, ifld
   character(len=8) char8, uname, vname, pre
   character(len=80) extrfile
   logical option,ex, vecflag
   
   if (trim(hfile%ftype)=='restart') then
      extrfile='extract.restart'
   else if (trim(hfile%ftype)=='nersc_daily') then
      extrfile='extract.daily'
   else if (trim(hfile%ftype)=='nersc_weekly') then
      extrfile='extract.weekly'
   else if (trim(hfile%ftype)=='archv') then
      extrfile='extract.archv'
   else if (trim(hfile%ftype)=='archm') then
      extrfile='extract.archm'
   else
      print *,'Unknown filetype '//trim(hfile%ftype)
      print *,'(fields_to_plot)'
      call exit(1)
   end if



   inquire(exist=ex,file=trim(extrfile)) 
   if (.not.ex) then
      print *,'Can not find extract file '//trim(extrfile)
      stop '(fields_to_plot)'
   end if


   open(21,file=trim(extrfile),status='old')
   read(21,*) ivers  
   if (ivers /= iversion) then
      close(21)
      print *,'Version of extract file should be:',iversion
      print *,'From 1 to 2 just add a F or T in column 3,5,7 of line 3'
      print *,'This makes it possible to plot solution on a sphere,'
      print *,'to calculate and plot rotated velocities, and to switch of'
      print *,'the index velocities. Use the line:'
      print *,'2 T T T     # 2-D, sphere, rotate, index velocities'
      print *
      print *,'From version 3 to 4, field names are b characters long, '
      print *,'starting on first column'
      stop 'wrong version of extract file'
   endif
   read(21,*) ! ignored
   read(21,*) ! ignored
   read(21,*) ! ignored
   read(21,*) ! ignored

   print *,'The following variables (and corresponding layers) will be printed:'
   nfld=0
   do i=1,1000
      read(21,*,end=888) char8,lay1,lay2,option  
      if (option) then
         nfld=nfld+1
         fld(nfld)%fextract=adjustl(char8)
         fld(nfld)%laya=lay1
         fld(nfld)%layb=lay2
         fld(nfld)%option=option
         print '(a5,i2,i2)',fld(nfld)%fextract,fld(nfld)%laya,fld(nfld)%layb
      endif
   enddo
888 close(21)

   ! Special case (and rules) for vectors -
   ! 1) No scalar variables should begin with u,v,taux or tauy
   ! 2) v(tauy)-component must come immediately after u(taux)-component
   ! TODO: improve this test
   fld(1:nfld)%vecflag=.false.
   do ifld=1,nfld
      if (fld(ifld)%fextract(1:1)=='u' .or. fld(ifld)%fextract(1:4)=='taux') then
      if (ifld==nfld) then
         print *,'For Vectors, the v-component must come immediately after the u component'
         print *,'(hint - dont place vars beginning with "u" at the end of the extract list)'
         print *,fld(ifld)
         print *,fld(ifld+1)
         stop
      else if (fld(ifld)%fextract(1:1)=='u' .and.  &
         (fld(ifld+1)%fextract(1:1)/='v' .or. fld(ifld)%fextract(2:8)/=fld(ifld+1)%fextract(2:8) )) then
         print *,'For Vectors, the v-component must come immediately after the u component - 1'
         print *,fld(ifld)
         print *,fld(ifld+1)
         stop
      else if (fld(ifld)%fextract(1:4)=='taux' .and. &
         ( fld(ifld)%fextract(5:8)/=fld(ifld+1)%fextract(5:8) .or. fld(ifld+1)%fextract(1:4)/='tauy' )) then
         print *,'For Vectors, the v-component must come immediately after the u component - 4'
         print *,fld(ifld)%fextract(1:4)
         print *,fld(ifld  )%fextract
         print *,fld(ifld+1)%fextract
         stop
      end if
      end if

      
      uname=fld(ifld  )%fextract
      vname=fld(ifld+1)%fextract ; 
      ! Flag indicating vector component
      vecflag= &
         (uname(1:1) == 'u'   .and.vname(1:1) == 'v') .or. &
         (uname(1:4) == 'taux'.and.vname(1:4) == 'tauy')  
      ! Now test that the rest of their names are equal
      if (vecflag) then
         if (uname(1:4) == 'taux') then
            vecflag=uname(5:8) == vname(5:8)
            pre='tau'//uname(5:8)
         elseif (uname(1:1) == 'u' .and. len_trim(uname)>1) then
            vecflag=uname(2:8) == vname(2:8)
            pre=uname(2:8)
         end if
      end if

      if (vecflag) then
         fld(ifld  )%vecflag = .true.
         fld(ifld+1)%vecflag = .true.
         fld(ifld  )%vecpost = pre
         fld(ifld+1)%vecpost = pre

         ! Switch off option for 2nd component, since it will be processed
         ! together with the 1st component
         fld(ifld+1)%option = .false.
         write(6,'(a)') 'Found vector pair: '//uname//vname
      end if

   end do


end subroutine fields_to_plot
end module m_fields_to_plot
