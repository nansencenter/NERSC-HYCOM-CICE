module m_ncvar_read

      integer, parameter, private :: lprt=6

contains

   subroutine ncvar_read(filename,varname,field,nx,ny,nrec,rec1,rec2,forcerecopt)
      use netcdf
      use mod_xc , only : xcstop,mnproc
      implicit none

      character(len=*),            intent(in)    :: varname,filename
      integer,                     intent(in)    :: nx,ny,nrec
      integer,                     intent(in)    :: rec1,rec2
      real, dimension(nx,ny,nrec), intent(inout) :: field
      integer, optional, intent(in)              :: forcerecopt


      integer :: i
      integer :: ncid,rdim,ngatts,varid,ndims,natts
      integer :: recdim,lvec
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: dimsizes(NF90_MAX_VAR_DIMS)
      integer :: alloc_size(NF90_MAX_VAR_DIMS)
      integer :: startvec(NF90_MAX_VAR_DIMS)
      integer :: mapvec  (NF90_MAX_VAR_DIMS)
      real :: add_offset, scale_factor
      real, allocatable :: tmp1(:), tmp2(:,:),tmp3(:,:,:)
      integer :: forcerec
      logical, parameter :: quiet=.true.

      forcerec=-100 ! dimensions IDs are from 0 (1?) and up, -100 should give no match
      if (present(forcerecopt)) forcerec=forcerecopt

      

      !print *,'ncvar_read entry'
      ! Error checks
      if (rec2-rec1> nrec ) then 
         if (mnproc==1) then
            write(lprt,'(a)') 'Attempted to extract more records than'// &
                              'place for in temporari field "field"  '
            write(lprt,'(a)') 'Filename is        : '//filename
            write(lprt,'(a)') 'Variable name is is: '//varname
         end if
         call xcstop('(ncvar_read)')
         stop '(ncvar_read)'
      else if (rec2<0.or.rec1<0) then
         if (mnproc==1) write(lprt,'(a)') 'Incorrect record indices rec1 rec2 ...'
         call xcstop('(ncvar_read)')
         stop '(ncvar_read)'
      else if (rec2-rec1+1<nrec) then
         if (mnproc==1) write(lprt,'(a)') 'Warning: "field" will not be filled up '
      end if


      ! Open NetCDF file
      call nf90_handle_err( &
           NF90_OPEN(filename,NF90_NOCLOBBER,ncid),&
           'Can not open '//trim(filename))


      ! Get unlimited (record) dimension and num of global atts
      call nf90_handle_err(&
           nf90_Inquire(ncid, unlimitedDimId=rdim,nattributes=ngatts),&
           'Error On inquiring on record dimension ')


      ! Output global attribute info if wanted
      !print *,'ncvar_dump_att 1 entry'
      call dump_atts(ncid,NF90_GLOBAL,ngatts,quiet=quiet)
      !print *,'ncvar_dump_att 1 exit'
               

      ! Inquire on variable -- get its id
      call nf90_handle_err( &
           nf90_inq_varid(ncid, varname, varid),   &
           filename//' does not contain variable '//varname)

     ! Inquire on variable -- dimensions and dimension ids
     call nf90_handle_err( &
          nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts),&
          'Error on inquiring variable '//varname)


     ! Error if ndim>3 (Max dimension of "field"
     if (ndims>3) then
        if (mnproc==1) &
        write(lprt,'(a)') 'Sorry -- this routine can only '// &
                          'handle 3 dimensonal Netcdf Vars'
        call xcstop('(ncvar_read)')
        stop '(ncvar_read)'
     end if
        

     ! Display variable atts -- return scale_factor and add_offset if
     ! necessary
      !print *,'ncvar_dump_att 2 entry'
      call dump_atts(ncid,varid,natts,scale_factor,add_offset,quiet)
      !print *,'ncvar_dump_att 2 exit'

     ! Inquire on variable dimensions -- sizes
     dimsizes=0
     alloc_size=0
     recdim=0
     lvec=1
     do i=1,ndims
        call nf90_handle_err( &
             nf90_Inquire_Dimension(ncid, dimids(i), len=dimsizes(i)), &
             ' Error on inquiring dimension for var '//varname)
       if (dimids(i)==rdim .or. i==forcerec)  then
          !recdim=rdim
          recdim=i
          lvec=lvec*nrec
          alloc_size(i)=nrec
          !print *,i,alloc_size(i),lvec,nrec
       else
          lvec=lvec*dimsizes(i)
          alloc_size(i)=dimsizes(i)
          !print *,i,alloc_size(i),lvec
       end if
     end do


     ! Check that field dimensions are consistent
     if (lvec /= nx*ny*nrec ) then
        if (mnproc==1) then
           write(lprt,'(a)') 'field dims are inconsistent'
           print *,'nx ny nrec ',nx,ny,nrec
           print *,'nx*ny*nrec ',nx*ny*nrec,lvec
        end if
        call xcstop('(ncvar_read)')
        stop '(ncvar_read)'
     else if (recdim>0) then
        !print *,'record dim : ',dimsizes(recdim)
        if  (rec1>dimsizes(recdim) .or. rec2>dimsizes(recdim)) then
           if (mnproc==1) write(lprt,'(a)') 'rec1 or rec2 exceed record dimension'
           call xcstop('(ncvar_read)')
           stop '(ncvar_read)'
        end if
     else if (recdim==0) then
        if (rec1/=1 .or. nrec>1 ) then
           if (mnproc==1) write(lprt,'(a)') 'Variable '//varname//' does not contain'//&
                 'a record id but nrec > 1'
           call xcstop('(ncvar_read)')
           stop '(ncvar_read)'
        end if
     end if


     ! Set "start" vector
     do i=1,ndims
          if (dimids(i)==rdim .or. i==forcerec)  then
            startvec(i)=rec1
         else
            startvec(i)=1
         end if
      end do



      ! Finally read field 
      !print *,startvec(1:ndims)
      !print *,dimsizes(1:ndims)
      !print *,alloc_size(1:ndims)
      !print *,nx,ny,nrec,rec1,rec2,lvec
      if (ndims==1) then
         allocate(tmp1(alloc_size(1)))
         call nf90_handle_err( &
              NF90_GET_VAR(ncid,varid,tmp1,&
              start=startvec(1:ndims)),                        &
              'Error in retrieving '//varname//' values')      
        field = reshape(tmp1,(/nx,ny,nrec/))
        !print *,field
        ! print *,scale_factor,add_offset,maxval(field),minval(field)
        field = field * scale_factor + add_offset
        deallocate(tmp1)

      else if (ndims==2) then
         allocate(tmp2(alloc_size(1),alloc_size(2)))
         call nf90_handle_err( &
              NF90_GET_VAR(ncid,varid,tmp2,&
              start=startvec(1:ndims)),                        &
              'Error in retrieving '//varname//' values')      
        field = reshape(tmp2,(/nx,ny,nrec/))
        !print *,scale_factor,add_offset,maxval(field),minval(field)
        field = field * scale_factor + add_offset
        deallocate(tmp2)

      else if (ndims==3) then
         allocate(tmp3(alloc_size(1),alloc_size(2),alloc_size(3)))
         call nf90_handle_err( &
              NF90_GET_VAR(ncid,varid,tmp3,&
              start=startvec(1:ndims)),                        &
              'Error in retrieving '//varname//' values')      
        field = reshape(tmp3,(/nx,ny,nrec/))
        !print *,scale_factor,add_offset,maxval(field),minval(field)
        field = field * scale_factor + add_offset
        deallocate(tmp3)
      end if


      ! Close NetCDF file
      call nf90_handle_err(&
           NF90_CLOSE(ncid),&
           'Error on closing '//filename)

      !print *,'ncvar_read exit'
   end subroutine ncvar_read







   subroutine dump_atts(ncid,varid,natts,scale_factor,add_offset,quiet)
      use netcdf
      implicit none
      integer,intent(in)  :: ncid,varid,natts
      real,intent(out),optional :: scale_factor,add_offset
      logical, intent(in), optional :: quiet

      logical :: silent
      integer :: i
      integer :: nctype,attlen
      character(len=80) :: attname
      character(len=3) :: cll
      character(len=200) :: attval_char
      real :: attval_real(200)

      silent = .true.
      if (present(quiet)) silent = quiet
      if(present(scale_factor)) scale_factor=1.
      if(present(add_offset)) add_offset=0.

      
      if (.not.silent) write(lprt,'(a)') 'Global NetCDF file info '
      do i=1,natts
         call nf90_handle_err(&
              NF90_INQ_ATTNAME(ncid,varid,i,attname), &
              'Error on inquiring attribute name')

         !print *,'attname = '//trim(attname)

         ! Skip history attribute (can be really long and causes (recoverable) error messages)
         if (trim(attname)=='history') cycle

         call nf90_handle_err(&
              NF90_INQUIRE_ATTRIBUTE(ncid,varid,attname,xtype=NCtype,len=attlen), &
              'Error on inquiring attribute '//attname)

         !print *,i,natts,ncid,varid,trim(attname),attlen,NCtype
         if (attlen<200) then
            if (NCtype==NF90_CHAR) then

               call nf90_handle_err(&
                    NF90_GET_ATT(ncid,varid,trim(attname),attval_char(1:attlen)),&
                    'Error obt att')

               if(.not.silent)&
                  write(lprt,'(a)') trim(attname)//':'//trim(attval_char(1:attlen))


            else if (NCtype==NF90_SHORT.or.NCtype==NF90_INT.or. &
                     NCtype==NF90_FLOAT.or.NCtype==NF90_DOUBLE ) then

               !print *,NCtype,attlen
               !print *,'int=',NF90_INT ,'  short=', NF90_SHORT ,'   float=', NF90_FLOAT &
               !,'   double=', NF90_DOUBLE


               call nf90_handle_err(&
                    NF90_GET_ATT(ncid,varid,trim(attname),attval_real(1:attlen)),&
                    'Error obt att')

               write(cll,'(i3.3)')attlen
               if(.not.silent)&
                  write(lprt,'(a,'//cll//'g14.7)') trim(attname)//':',attval_real(1:attlen)
            else
               if(.not.silent)&
                  write(lprt,'(a,i4,a)') 'Variable type ',NCtype,' is unknown'
            end if

            if (attname=='scale_factor') then
               if (present(scale_factor)) &
                  call nf90_handle_err(&
                       NF90_GET_ATT(ncid,varid,attname,scale_factor), &
                       'Error obt att')
            else if (attname=='add_offset') then
               if (present(add_offset))  &
               call nf90_handle_err(&
                    NF90_GET_ATT(ncid,varid,attname,add_offset),&
                 'Error obt att')
            end if

         else
            !print *,attlen
            if(.not.silent)&
            print *, 'Length of attribute '//trim(attname)//' > 200, I will not display it ...'
         end if
      end do
   end subroutine




   subroutine nf90_handle_err(errcode,info)
      use mod_xc , only: xcstop, mnproc
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: info

      if (errcode/=NF90_NOERR) then
         if (mnproc==1) then
            write(lprt,'(a)') NF90_STRERROR(errcode)
            write(lprt,'(a)') info
         end if
         stop '(ncvar_read)'
         call xcstop('(ncvar_read)')
      end if
   end subroutine



end module m_ncvar_read














      
         

      

