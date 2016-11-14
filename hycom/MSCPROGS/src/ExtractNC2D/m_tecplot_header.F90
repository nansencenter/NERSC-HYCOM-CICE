module m_tecplot_header
contains
subroutine tecplot_header(normal,sphere,rotate,nz,fld,nfld,runcode)
   use mod_types
   implicit none
   integer, intent(in)          :: nz
   integer, intent(in)          :: nfld
   type(fields),intent(in)      :: fld(nfld)
   logical, intent(in)          :: normal
   logical, intent(in)          :: sphere
   logical, intent(in)          :: rotate
   character(len=3), intent(in) :: runcode

   character(len=20)  string
   character(len=2)  clay
   integer i,k

   !print *,'Tecplot header:'

   open(10,file='head.txt')
      write(10,'(a,a3,a)') 'TITLE = "',runcode,'"'
      write(10,'(a)',advance='no')'VARIABLES = "I-index" "J-index"'
      write(10,'(a)',advance='no')'"Lon" "Lat" "Depth" '
      if (sphere) then
         write(10,'(a)',advance='no')' "X-new" "Y-new" "Z-new" '
      endif

! printing all fields (or variable) names in tecplot header.
! Velocity names only plotted if normal is true!!!!

      do k=0,nz
      do i=1,nfld
         if ((fld(i)%laya <= k) .and. (k <= fld(i)%layb).and. &
             fld(i)%option .and. (.not. fld(i)%vecflag)) then
            string(1:1)='"'
            string(2:9)=fld(i)%fextract 
            string(10:10)='('
            write(string(11:12),'(i2.2)')k
            string(13:15)=')" '

            if (fld(i)%fextract(1:1) == 'U'.or.fld(i)%fextract(1:1) == 'V' .or.  &
                fld(i)%fextract(1:1) == 'u'.or.fld(i)%fextract(1:1) == 'v') then
               if (normal) then
                  write(10,'(a20)',advance='no')string(1:15)
                  !print  *,string(1:12)
               endif
            else
               write(10,'(a20)',advance='no')string(1:15)
               !print  *,string(1:12)
            endif

         endif
      enddo
      enddo


      if (rotate .or. sphere) then
         do k=1,nz
            do i=1,nfld-1
               if (((fld(i)%laya <= k) .and. (k <= fld(i)%layb)) .and.  &
                    (fld(i)%option .and. fld(i)%vecflag)) then

                  if (normal) then
                     write(clay,'(i2.2)')k
                     write(10,'(a18)',advance='no') '"'//fld(i)%fextract//clay//'      "'
                     write(10,'(a18)',advance='no') '"'//fld(i)%fextract//clay//'      "'
                  end if

                  if (rotate) then

                     string(2:6)='UROT_'
                     string(7:13)=fld(i)%fextract(2:8)
                     string(14:14)='('
                     write(string(15:16),'(i2.2)')k
                     string(17:18)=')"'

                     write(10,'(a18)',advance='no')string(1:18)
                     !write(*,'(a18)',advance='no')string(1:18)

                     string(2:2)='V'
                     write(10,'(a18)',advance='no')string(1:18)
                     !write(*,'(a18)',advance='no')string(1:18)
                  endif

                  if (sphere) then
                     string(2:6)='USPH_'
                     string(7:13)=fld(i)%fextract(2:8)
                     string(14:14)='('
                     write(string(15:16),'(i2.2)')k
                     string(17:18)=')"'

                     write(10,'(a18)',advance='no')string(1:18)

                     string(2:2)='V'
                     write(10,'(a18)',advance='no')string(1:18)

                     string(2:2)='W'
                     write(10,'(a18)',advance='no')string(1:18)
                  endif
                  !write(*,*)
                  write(10,*)
               endif
            enddo
         enddo
         write(10,'(a1)')' '
      endif

   close(10)
   !print '(a)','#############################################'

end subroutine tecplot_header
end module m_tecplot_header
