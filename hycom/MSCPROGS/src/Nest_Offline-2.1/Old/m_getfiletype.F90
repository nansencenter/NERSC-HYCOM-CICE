module m_getfiletype
contains

      subroutine getfiletype(filename,fbase,hycfile)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=80), intent(out) :: fbase
      integer, intent(out) :: hycfile
      integer :: findhdr,findab,finddaily,findweek,findrst

      ! Check for type ...
      findhdr  =index(filename,'.hdr')
      findab=max(index(filename,'.a'),index(filename,'.b'))
      if (.not. findab>0 .and. .not. findhdr>0) then ! old pak type

         print *,'No .ab or .hdr files'
         hycfile=-1

      else

         ! We have a .ab-file. Now figure out what type.
         findrst  =index(filename,'restart')
         finddaily=index(filename,'DAILY')
         findweek =index(filename,'AVE')
         findhdr  =index(filename,'.hdr')
         print *,findhdr
         if (findrst==4) then
            hycfile=1 ! 1 for restart files
            fbase=filename(1:findab-1)
         elseif (finddaily==4) then
            hycfile=2 ! 1 for daily files
            fbase=filename(1:findab-1)
         elseif (findweek==4) then
            hycfile=3 ! 1 for weekly files
            fbase=filename(1:findab-1)
         elseif (findhdr>0) then
            hycfile=4 ! 1 for weekly files
            fbase=filename(1:findhdr-1)
         else
            print *,'Can not deduce file type from  file name'
            hycfile=-1
         end if
      end if

      end subroutine getfiletype

end module m_getfiletype
