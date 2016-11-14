module m_get_hyc_infiles
contains
subroutine get_hyc_infiles(hyc_files,numfiles,maxfiles)
   implicit none 

   integer , intent(out) :: maxfiles, numfiles
   character(len=80), dimension(:),pointer :: hyc_files


   character(len=*), parameter  :: infile_hyc    = 'infiles.hyc'
   integer         , parameter  :: l_maxfiles=1000

   logical :: ex
   integer :: ios, counter

   print *
   print '(a)','Reading infiles info:'

   ! Open infile for program -- files to process
   inquire(exist=ex,file=infile_hyc)
   if (.not. ex) then
      print '(a)','You must specify files to process in'
      print '(a)','the file  '//infile_hyc
      stop '(get_hyc_infiles)'
   endif


   ! Cycle through infile to get files to process
   ios=0
   counter=1
   maxfiles=l_maxfiles
   allocate(hyc_files(maxfiles))
   open(10,file=infile_hyc,form='formatted')
   do while (ios==0 .and. counter<maxfiles)
      read(10,*,iostat=ios) hyc_files(counter)
      if (ios==0) then
         print '(a)','INPUT FILE: '//trim(hyc_files(counter))
         counter=counter+1
      end if
   end do
   numfiles=counter-1
   if (counter==maxfiles) then
      print *,'WARNING: max num of files reached !'
      stop '(get_hyc_files)'
   end if
   print '(a,i4,a)','Read ',numfiles,' files'
   close(10)
end subroutine get_hyc_infiles
end module m_get_hyc_infiles
