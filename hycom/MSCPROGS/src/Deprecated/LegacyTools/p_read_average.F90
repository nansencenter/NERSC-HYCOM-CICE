!KAL - routine for reading old "weekly average" file types - no other processing is done here
!KAL - grid size is hardcoded in mod_dimensions
program read_ave
   use mod_average
   implicit none
   logical ex
   character(len=80) :: tmpchar
   integer *4 , external :: iargc


   if (iargc()==1) then
      call getarg(1,tmpchar)
   else
      print *,'Supply file name at command line'
      stop '(read_average)'
   end if

   print *,'reading:', trim(tmpchar)
   inquire(file=trim(tmpchar),exist=ex)
   if (.not.ex) then
      print *,'File does not exist'
      stop '(read_average)'
   end if
   open(10,file=trim(tmpchar),status='unknown',form='unformatted')
   read(10)ave_week
   close(10)
   print *,'counter= ',ave_week%counter

   ! Example layer variables t s etc are layer weighted, d is layer thickness
   print *,'min temp layer 1',minval(ave_week%t(:,:,1)/ave_week%d(:,:,1),mask=ave_week%d(:,:,1)>1   )
   print *,'max temp layer 1',maxval(ave_week%t(:,:,1)/ave_week%d(:,:,1),mask=ave_week%d(:,:,1)>1   )
   print *,'max dp   layer 1',maxval(ave_week%d(:,:,1)/ave_week%counter ,mask=ave_week%d(:,:,1)>1)


   ! Example 2D variables
   print *,'min ssh ',minval(ave_week%ssh(:,:)/ave_week%counter,mask=ave_week%d(:,:,1)>1)
   print *,'max ssh ',maxval(ave_week%ssh(:,:)/ave_week%counter,mask=ave_week%d(:,:,1)>1)

end program

