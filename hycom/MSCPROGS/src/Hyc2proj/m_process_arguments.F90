module m_process_arguments
contains

subroutine process_arguments(cmethod,listmode,filestart)
implicit none
character(len=*),  intent(out) :: cmethod
logical,  intent(out) :: listmode
integer,  intent(out) :: filestart
character(len=80) :: tmparg
logical :: lerr
integer :: i
#if defined (IARGC)
integer*4, external :: iargc
#endif

   ! Test for arguments
   lerr=.false.
   cmethod=''
   filestart=1
   if (iargc()>0) then
      lerr=.false.

      call getarg(1,tmparg)

      ! Denotes vertical interpolation method
      if (trim(tmparg)=='--vertint') then
         if (iargc()<3) then
            lerr=.true.
         else
            call getarg(2,cmethod)
            if ( trim(cmethod)=='spline' .or.  trim(cmethod)=='staircase' .or.  &
                 trim(cmethod)=='linear' ) then
                 print *,'Vertical interpolation:'//trim(cmethod)
                 filestart=3
            else
               lerr=.true.
            end if
         end if

      ! Denotes other options (presently just --help , which triggers error
      ! anyway...)
      else if (tmparg(1:2)=='--') then
         lerr=.true.
      end if
   else
      lerr=.true.
   end if

   if (lerr) then 
      print *,'Routine hyc2proj is used to interpolate from HYCOM model fields'
      print *,'onto a grid projection and at specified depth levels. To use this ' 
      print *,'routine several input files are needed. You need a file to specify' 
      print *,'the projection, "proj.in". You need a file to specify depth levels,' 
      print *,'"depthlevels.in", and finally you need a file to specify fields to '
      print *,'extract (extract.restart, extract.daily, extract.weekly,'
      print *,'extract.archv, extract.archv_wav). Files to process are given on the'
      print *,'command line.'
      print *,'Vertical interpolation method can also be given on command line,'
      print *,'defaul method is spline.'
      print *
      print *,'Options are :'
      print *,' hyc2proj [--vertint value] file1 file2 ...'
      print *,'   where value is either spline, linear or staircase'
      print *,'   staircase uses layer value at specified depths          (fast)'
      print *,'   linear    uses linear weighting between layer midpoints (fast)'
      print *,'   spline fits a spline to the depth profile               (slow)'
      print *,'(process_arguments)'
      call exit(1)
   end if

   if (trim(cmethod)=='') then
      cmethod='spline'
   end if

   end subroutine
end module
