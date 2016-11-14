module m_read_blkdat
contains

   !Reads different info from blkdat.input 

   subroutine read_gridindex(idm,jdm,kdm)
      implicit none

      ! Dummy variables
      integer, intent(out) :: idm,jdm,kdm
      integer    dump1,dump2 
      integer iversn,iexpt
      character(len=80) :: file1
      logical :: ex
      character(len=80) :: ctitle(4)


      ! Check for existence of blkdat.input file
      file1='blkdat.input'
      inquire(file=trim(file1),exist=ex)
   if (.not.ex) then
      print *,'File blkdat.input does not exist'
      stop
   end if
 
    

open(10,file=file1,form='formatted',status='old')
 read(10,116) ctitle,iversn,iexpt,idm,jdm,dump1,dump2,kdm

close (10)
!116 format (a80/a80/a80/a80/ &
!      i4,4x,'''iversn'' = hycom version number x10'/  &
!      i4,4x,'''iexpt '' = experiment number x10'/  &
!      i4,4x,'''idm   '' = longitudinal array size'/  &
!      i4,4x,'''jdm   '' = latitudinal  array size'/  &
!      i4,4x,'''itest '' = Year of integration start '/  &
!      i4,4x,'''jtest '' = Day of integration start'/  &
!      i4,4x,'''kdm   '' = Vertical     array size')  
116  format (a80/a80/a80/a80/ &
       i4,4x,35x/ &
       i4,4x,32x/  &
       i4,4x,34x/  &
       i4,4x,34x/  &
       i4,4x,37x/  &
       i4,4x,35x/  &
       i4,4x,34x)  

   end subroutine read_gridindex

   subroutine read_thbase(thbase)
      implicit none

      ! Dummy variables
      real, intent(out) :: thbase
      integer    dump1,dump2 
      integer iversn,iexpt
      character(len=80) :: file1
      logical :: ex
      character(len=80) :: ctitle(22)


      ! Check for existence of blkdat.input file
      file1='blkdat.input'
      inquire(file=trim(file1),exist=ex)
   if (.not.ex) then
      print *,'File blkdat.input does not exist'
      stop
   end if
 
    

open(10,file=file1,form='formatted',status='old')
 read(10,116) ctitle,thbase
print*,'thbase inside read_thbase',thbase
close (10)
116  format (a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/ &
             a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/ &
       2x,f3.1,3x,'''thbase'' = reference density (sigma units)') 

   end subroutine read_thbase


end module m_read_blkdat
