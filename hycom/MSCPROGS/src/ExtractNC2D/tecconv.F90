program tecconv
   use mod_xc
   use mod_za
   use mod_grid
   use mod_types
   use mod_hycomfile_io, only: hycomfile, initHF, getfiletype, vDim
   use netcdf
   use m_handle_err
   use m_fields_to_plot
   use m_tecplot_header
   use m_tecplot_zoneinfo
   use m_tecplot_dump
   use m_tecplot_dump_rot
   implicit none
   type(fields)       fld(1000)
   character(len=3)   rungen
   integer nfld, i

   logical sphere,rotate,normal,ex
   character(len=8) dpfield,dptype
   integer :: rdimid,jdimid,idimid,ncid,nrec
   integer*4, external :: iargc
   character(len=20) :: tmparg
   character(len=24) :: plot_time

   character(len=80) :: filename, ftype
   integer :: kdm
   real    :: rtime
   real, dimension(:,:), allocatable :: plon2, plat2
   type(hycomfile) :: hfile

   ! KAL - not all these are implemented yet
   logical :: ltecplot     ,   &   ! Input flag -- produce tecplot output
              skipheader           ! TODO: Input flag -- skip header in .b files

   ltecplot=.false.
   skipheader=.false.


   if (iargc()==1) then
      call getarg(1,filename) 
   elseif (iargc()>=2) then
      do i=1,iargc()-1
         call getarg(i,tmparg) 

         if (trim(tmparg) == '-tec') then
            ltecplot=.true.
         elseif (trim(tmparg) == '-weekly') then
            ftype='nersc_weekly'
         elseif (trim(tmparg) == '-daily') then
            ftype='nersc_daily'
         elseif (trim(tmparg) == '-restart') then
            ftype='restart'
         elseif (trim(tmparg) == '-skipheader') then
            skipheader=.true.
            ! ingored for now
         end if
      end do
      call getarg(iargc(),filename) 
   else
      print *,'file name must be supplied '
      stop '(tecconv)'
   end if
   rungen = filename(1:3)

   ! Init io modules
   call xcspmd()
   call zaiost()

   ! What file is this? (daily? weekly? restart? pak?)
   ! TODO: force filetype option goes here
   ftype=getfiletype(trim(filename))

   ! Inits file type
   call initHF(hfile,trim(filename),trim(ftype))
   rtime=hfile%fyear
   kdm=vDim(hfile)

   ! Get model grid
   call get_grid()

! determine number of and which fields to be plotted (stored in fld), 
! the number of fields (nfld) and options for grid projections.
   call fields_to_plot(sphere,rotate,normal,fld,nfld,hfile,kdm)
   

!  Write header for tecplot file
   if (ltecplot) then
      call tecplot_header(normal,sphere,rotate,kdm,fld(1:nfld),nfld,rungen)
   end if


! Write zone header and position of x-, y-, and z-points
   write(plot_time,'("d",i4.4," d",i3.3," h",i2.2)') hfile%iyear, hfile%iday, hfile%ihour
   !plot_time='d'//rtd%cyy//' m'//rtd%cmm//' d'//rtd%cdm//' h'//rtd%chh
   allocate(plon2(0:idm+1,0:jdm+1))
   allocate(plat2(0:idm+1,0:jdm+1))
   call lon_lat_extrap(plat,plon,plat2,plon2,idm,jdm)
   call tecplot_zoneinfo(idm,jdm,kdm,plot_time,qlon(1:idm+1,1:jdm+1),qlat(1:idm+1,1:jdm+1), &
                         scuy,scvx,depths,sphere,ltecplot)

! Get record number for Netcdf files
   if (NF90_OPEN('tmp1.nc',NF90_WRITE,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(restart2netcdf)'
   end if
   call handle_err (nf90_inq_dimid(ncid, 'idim', idimid))
   call handle_err (nf90_inq_dimid(ncid, 'jdim', jdimid))
   call handle_err (nf90_inq_dimid(ncid, 'rdim', rdimid))
   call handle_err (nf90_Inquire_Dimension(ncid, rdimid, len=nrec))
   nrec=nrec+1
   call handle_err(NF90_CLOSE(ncid))


! Find, extract and dump all scalar variables 
   call tecplot_dump(depths,idm,jdm,kdm,fld(1:nfld),nfld,normal,nrec,ltecplot,hfile)

! Find, extract and dump velocities rotated to lat-lon and/or projected on a sphere'
   if (rotate .or. sphere .or. normal) then
      call tecplot_dump_rot(depths,idm,jdm,kdm,fld,nfld,plat2,plon2,normal,rotate,sphere,nrec,ltecplot,hfile)
   endif


end program tecconv
