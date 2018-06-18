module m_parse_blkdat
!   private :: blkini, blkinr, blkinvoid,  &
   private :: blkinvoid,  &
      parse_blkdatmain, parse_blkdatreal, &
      parse_blkdatint

   interface parse_blkdat
      module procedure parse_blkdatreal, &
                       parse_blkdatint
   end interface
contains


   

     
      subroutine parse_blkdatreal(cvar,realvar,blkfilein,imatch)
      implicit none
      character(len=6), intent(in)  :: cvar
      real   ,          intent(out) :: realvar
      character(len=*), intent(in), optional :: blkfilein
      integer         , intent(in), optional :: imatch
      integer :: intvar, imatch2
      character(len=80) :: blkfile2

      if (present(blkfilein)) then
         blkfile2=trim(blkfilein)
      else
         blkfile2='blkdat.input'
      end if

      if (present(imatch)) then
         imatch2=imatch
      else
         imatch2=1
      end if
      call parse_blkdatmain(cvar,'real',realvar,intvar,trim(blkfile2),imatch2)
      end subroutine

     
      subroutine parse_blkdatint(cvar,intvar,blkfilein,imatch)
      implicit none
      character(len=6), intent(in)  :: cvar
      integer,          intent(out) :: intvar
      character(len=*), intent(in), optional :: blkfilein
      integer         , intent(in), optional :: imatch
      integer :: imatch2
      character(len=80) :: blkfile2
      real :: realvar

      if (present(blkfilein)) then
         blkfile2=trim(blkfilein)
      else
         blkfile2='blkdat.input'
      end if

      if (present(imatch)) then
         imatch2=imatch
      else
         imatch2=1
      end if
      call parse_blkdatmain(cvar,'integer',realvar,intvar,trim(blkfile2),imatch2)
      end subroutine


   
      subroutine parse_blkdatmain(cvar,vtype,realvar,intvar,blkfilein,imatch)
      implicit none
      character(len=6), intent(in)  :: cvar
      character(len=*), intent(in)  :: vtype
      integer,          intent(out) :: intvar
      real   ,          intent(out) :: realvar
      character(len=*), intent(in), optional :: blkfilein
      integer         , intent(in), optional :: imatch

      character(len=80) :: blkfile

      logical :: found,ex
      integer :: nmatch,imatch2

      if (present(blkfilein)) then
         blkfile=blkfilein
      else
         blkfile='blkdat.input'
      end if
      if (present(imatch)) then
         imatch2=imatch
      else
         imatch2=1
      end if



      inquire(exist=ex,file=trim(blkfile))

      nmatch=0
      if (ex) then
         open(99,file=trim(blkfile))


         !! Skip header
         read(99,*)
         read(99,*)
         read(99,*)
         read(99,*)

         found=.false.

         do while (.not.found)
            found = blkinvoid(cvar)

            if (found) then
               nmatch=nmatch+1
               !print *,found,nmatch,imatch2
               found=found.and.nmatch==imatch2
               !print *,found
            end if

         end do

         ! if found, read..
         if (found) then
            backspace(99)
            if (trim(vtype)=='integer') then
               call blkini(intvar,cvar)
            elseif (trim(vtype)=='real') then
               call blkinr(realvar,cvar,'(a6," =",f10.4," m")')
            else
               print *,'Dont know how to handle variable type '//trim(vtype)
               stop '(parse_blkdat)'
            end if
         else
            print *,'Cant find varable'
            stop '(parse_blkdat)'
         end if

         close(99)
      else
         print *,'Cant find '//trim(blkfile) 
         stop '(parse_blkdat)'
      end if
      end subroutine parse_blkdatmain




      subroutine blkinr(rvar,cvar,cfmt)
      !use mod_xc  ! HYCOM communication interface
      implicit none
      real      rvar
      character cvar*6,cfmt*(*)
!     read in one real value
      character*6 cvarin

      read(99,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
      call flush(6)

      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin, &
                            ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop '(blkinr)'
      endif
      return
      end subroutine

      subroutine blkini(ivar,cvar)
      implicit none
      integer     ivar
      character*6 cvar
!     read in one integer value
      character*6 cvarin
 
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
      call flush(6)
 
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin, &
                            ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop '(blkini)'
      endif
      return
 6000 format(a6,' =',i6)
      end subroutine


      subroutine blkinl(lvar,cvar)
      implicit none
      logical     lvar
      character*6 cvar
!     read in one logical value
      character*6 cvarin
      integer     ivar

      read(99,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
      call flush(6)

      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinl - input ',cvarin, &
                            ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop('(blkinl)')
      endif
      return
 6000 format(a6,' =',l10)
      end subroutine





      logical function blkinvoid(cvar)
      !use mod_xc  ! HYCOM communication interface
      implicit none
      real      rvar
      character cvar*6
!     read in one real value
      character*6 cvarin

      read(99,*) rvar,cvarin
      !print *,rvar,cvarin,cvar
      blkinvoid=trim(cvar)==trim(cvarin)
      end function


end module m_parse_blkdat
