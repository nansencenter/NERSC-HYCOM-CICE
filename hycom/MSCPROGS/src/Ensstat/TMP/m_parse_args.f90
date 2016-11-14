module m_parse_args
integer, parameter :: maxfiles=200
integer,save :: nfiles
character(len=8), save :: cfld1, cfld2,mode
integer         , save :: k1, k2
integer         , save :: xpoint, ypoint
character(len=80),save :: filelist(maxfiles)

contains
subroutine parse_args(idm,jdm,kdm)
implicit none
integer, intent(in) :: idm, jdm, kdm
integer*4, external :: iargc
character(len=80) :: tmparg
integer :: ifile
logical :: error,ex


   ! Command line argument - filename 
   error =.false.
   if (iargc()<1) then
      error=.true.
      goto 200
   end if



   call getarg(1,mode)
   if (trim(mode)=="--field") then

      if (iargc()>=6) then
         print *,"hei",mode
         call getarg(2,cfld1)
         call getarg(3,tmparg) ; read(tmparg,*) k1
         call getarg(4,cfld2)
         call getarg(5,tmparg) ; read(tmparg,*) k2
         if ( k1<0 .or.  k1 > kdm .or. k2<0 .or. k2>kdm) then
            print *,'vertical coordinate out of range'
            print *,'k1=',k1
            print *,'k2=',k2
            print *,'kdm=',kdm
            call exit(1)
         end if
         nfiles=iargc()-5
         do ifile=1,nfiles
            call getarg(ifile+5,filelist(ifile))
            inquire(exist=ex,file=filelist(ifile))
            if (.not. ex) then
               print *,'File '//trim(filelist(ifile))//' does not exist'
               call exit(1)
            end if
         end do
      else 
         error=.true.
         goto 200
      end if

   else if (trim(mode)=="--point") then
      if (iargc()>=8) then
         call getarg(2,cfld1)
         call getarg(3,tmparg) ; read(tmparg,*) k1
         call getarg(4,cfld2)
         call getarg(5,tmparg) ; read(tmparg,*) k2
         if ( k1<0 .or.  k1 > kdm .or. k2<0 .or. k2>kdm) then
            print *,'vertical coordinate out of range'
            print *,'k1=',k1
            print *,'k2=',k2
            print *,'kdm=',kdm
            call exit(1)
         end if
         call getarg(6,tmparg) ; read(tmparg,*) xpoint
         call getarg(7,tmparg) ; read(tmparg,*) ypoint
         if ( xpoint<1 .or.  xpoint > idm .or. ypoint<1 .or. ypoint>jdm) then
            print *,'point coordinate out of range'
            print *,'i  ,j   ',idm,jdm
            print *,'idm,jdm=',xpoint,ypoint
            call exit(1)
         end if
         nfiles=iargc()-7
         do ifile=1,nfiles
            call getarg(ifile+7,filelist(ifile))
            if (.not. ex) then
               print *,'File '//trim(filelist(ifile))//' does not exist'
               call exit(1)
            end if
         end do
      else 
         error=.true.
         goto 200
      end if
   else
      error=.true.
      goto 200
   end if


200 continue
   if (error) then
      print *,'Usage:'
      print *,'1) '
      print *,'  ensstat --field fld1 k1 fld2 k2 restartfiles'
      print *,'2) '
      print *,'  ensstat --point fld1 k1 fld2 k2 cpoint ypoint restartfiles'
      stop '(ensstat)'
   end if
end subroutine
end module
