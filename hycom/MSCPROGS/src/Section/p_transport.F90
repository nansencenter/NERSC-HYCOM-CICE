program transport_calc
   use mod_xc
   use mod_za
   use mod_grid
   use mod_sections
   use mod_transport
   use mod_hycomfile_io
   implicit none
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   character(len=80) :: fnamein, tmparg
   character(len=20) :: ftype
   type(hycomfile) :: hfile
   integer :: kdm,i
   logical :: append

   append=.false.
   if (iargc()>=1) then
      call getarg(1,fnamein) ! Always filename
      do i=2,iargc()
         call getarg(i,tmparg)
         if (trim(tmparg)=='-append') append=.true.
      end do
   else
      print *,'No input given '
      print *,'(section_transport)'
      call exit(1)
   end if


   ! What file is this? (daily? weekly? restart? pak?)
   ftype=getfiletype(trim(fnamein))

   ! Initialize hycom file
   call initHF(hfile,fnamein,trim(ftype))

   ! For safety
   kdm=vDim(hfile)
   if (kdm <1 .or. kdm > 100 )  then
      print *,'kdm is ',kdm
      print *,'This is probably a bug ...'
      call exit(1)
   end if

   ! Initialize IO for .ab files
   CALL XCSPMD()  ! -- Requires "regional.grid.b" to be present
   CALL ZAIOST()

   ! Get model grid & Depths 
   call get_grid()

   ! Read nodes along section -- this assumes "section_intersect" is already run
   call read_section_nodes()

   ! Calculate transport -- dump to files
   call transport(hfile,appendfile=append)

   write(6,'(a)') 'section_transport finished -- OK'
end program transport_calc
