!      Generate SSH fields directly from the old restart file.
!      you have to also have a BOT file 
!      by Fanf 10/06/2004
!      
!      Changed by "knutali" 20 Oct 2006 -- adapted to standard 
!      hycom restart files
!
!      KAL -- March 2007, calculates SLA as well as SSH, requires meanssh  file to be    present



    program p_ssh_from_state
       use mod_za
       use mod_xc
       use mod_read_rstab
       use mod_year_info
       use mod_sigma
       use m_parse_blkdat
       use netcdf
       implicit none

       type(year_info) :: rtdump
       
       character(len=9) :: rident
       character(len=80):: filename, filenameSSH,depthfile, fbase, matfile
       character(len=80) :: tmpchar
       character(len=80) :: tmp
       character(len=3)  :: rungen,cmem
       character(len=6)  :: cvarin
       character(len=11) :: tag11
       integer*4, external :: iargc
    
       integer :: i, j, k, reclA, ios, reclSSH, ios2
       logical :: ex
       integer :: kdm
       integer :: idm_dummy, jdm_dummy, indxa,indxb

      real, allocatable, dimension(:,:) ::   &
         temp, saln, th3d, dp, psikk, thkk, oneta, pbavg, ssh,  &
         depths,modlon, modlat, meanssh, sla
      real*8, dimension(:,:), allocatable :: iofld
      real*4, dimension(:,:), allocatable :: iofld4
      integer, allocatable, dimension(:,:) ::  imask
      real, allocatable, dimension(:,:,:) ::   &
         p, montg, thstar

      real, parameter :: onem  = 9806.
      !real, parameter :: thref = 1e-3
      !real, parameter :: radian= 57.2957795
      !real, parameter :: pi    = 3.1415926536
      real, parameter :: undef = -1e14
      real    :: thbase,rvar
      integer :: kapflg, thflag
      logical :: tbaric
      real    :: hmin, hmax
      integer :: nrmem, imem
      integer :: idummy
      real    :: rdummy

      real :: tmpth3d, tmpkapf
      character(len=80) :: ncfile
      integer :: ncid,ierr,idmid,jdmid,var_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined (MATLAB)
! KAL -- For dumping to matlab file
#include </export/fimm/local/Matlab-R14sp3/extern/include/fintrf.h>
      integer :: iret
      MWPOINTER :: mxCreateNumericMatrix, mxGetPr, mxClassIDFromClassName, matopen,  &
         mxCreateDoubleMatrix, matPutVariableAsGlobal, mp, pa1
      integer matputvariable, matclose
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       !include 'stmt_funcs.h'
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Specify run id and time on command line
      if (iargc()/=1 .and. iargc()/=2) then
         print *,'Usage : '
         print *,'ssh_from_state <filename> [member]'
         print *
         print *,'Where: '
         print *,'   <filename> is restart file name'
         print *,'   <member>   is ensemble member (optional, default=1)'
         print *,' member is needed to put sla/ssh into correct file record'
         print *, '(ssh_from_state)'
         call exit(1)
      else if (iargc() >= 1 ) then
         call getarg(1,filename) ; rungen = trim(filename)

         indxa=index(filename,'.a')-1
         indxb=index(filename,'.b')-1
         if (indxa<=0 .and. indxb<=0) then
            print *,'Can not get file base for file name'
            stop
         end if
         fbase=filename(1:max(indxa,indxb))

         if (iargc()==2) then
            call getarg(2,tmp) ;
            read(tmp,'(i3)') imem
         else
            imem=1
         end if
      end if
      print *,iargc()
      print *,trim(fbase)



      !print *,'Warning: tbaric is hardcoded'
      !print *,'Warning: thbase is hardcoded'
      !tbaric=.false.
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff for model
       
      call xcspmd
      call zaiost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get # layers from model - from header

      call rst_read_header(fbase,rtdump,nrmem,idm_dummy,jdm_dummy,kdm)
      if (idm/=idm_dummy .or. jdm/=jdm_dummy ) then
         print *, 'Mismatch between restart grid size and depths grid size';
         stop '(ssh_from_state)'
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read bathymetry and grid from regional.depths file

      allocate(imask  (idm,jdm))
      allocate(depths (idm,jdm))
      allocate(modlon (idm,jdm))
      allocate(modlat (idm,jdm))
      call zaiopf('regional.depth.a','old',99)
      call zaiord(depths, imask,.false., hmin,hmax,  99)
      call zaiocl(99)
      call zaiopf('regional.grid.a','old',99)
      call zaiord(modlon, imask,.false., hmin,hmax,  99)
      call zaiord(modlat, imask,.false., hmin,hmax,  99)
      call zaiocl(99)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! Use the parse_blkdat routine to get data from blkdat.input
      call parse_blkdat('kapref','integer',rdummy,kapflg)
      call parse_blkdat('thflag','integer',rdummy,thflag)
      call parse_blkdat('thbase','real',thbase,idummy)



      !inquire(file='blkdat.input',exist=ex)
      !if (.not.ex) then
      !   write(lp,*) 'ERROR: blkdat.input does not exist'
      !   stop '(ssh_from_state:blkdat.input parse)'
      !endif
      !open(10,file='blkdat.input')
      !cvarin=''
      !read(10,*)
      !read(10,*)
      !read(10,*)
      !read(10,*)

      !do while (cvarin/='kapflg')
      !   read(10,*) rvar,cvarin
      !end do
      !if (cvarin/='kapflg') then
      !   write(lp,*) 'ERROR: Could not get jdm'
      !   stop '(ssh_from_state:blkdat.input parse)'
      !else
      !   kapflg=nint(rvar)
      !end if

      !do while (cvarin/='thflag')
      !   read(10,*) rvar,cvarin
      !end do
      !if (cvarin/='thflag') then
      !   write(lp,*) 'ERROR: Could not get thbase'
      !   stop '(ssh_from_state:blkdat.input parse)'
      !else
      !   thflag=nint(rvar)
      !end if

      !do while (cvarin/='thbase')
      !   read(10,*) thbase,cvarin
      !end do
      !if (cvarin/='thbase') then
      !   write(lp,*) 'ERROR: Could not get thbase'
      !   stop '(ssh_from_state:blkdat.input parse)'
      !end if

      tbaric=kapflg==-1
      print *,'blkdat.input -- kapflg:',kapflg
      print *,'blkdat.input -- thflag:',thflag
      print *,'blkdat.input -- thbase:',thbase
      print *,'blkdat.input -- tbaric:',tbaric


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate fields

      allocate(iofld(idm,jdm))
      allocate(iofld4(idm,jdm))
      allocate(temp  (idm,jdm))
      allocate(saln  (idm,jdm))
      allocate(th3d  (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(psikk (idm,jdm))
      allocate(thkk  (idm,jdm))
      allocate(oneta (idm,jdm))
      allocate(pbavg (idm,jdm))
      allocate(ssh   (idm,jdm))
      allocate(sla   (idm,jdm))
      allocate(meanssh(idm,jdm))
      allocate(p     (idm,jdm,kdm+1))
      allocate(montg (idm,jdm,kdm))
      allocate(thstar(idm,jdm,kdm))




      ! Get restart fields psikk and thkk
      call read_rstfield2d(fbase,'psikk   ',psikk,idm,jdm,0,undef)
      call read_rstfield2d(fbase,'thkk    ',thkk ,idm,jdm,0,undef)
      call read_rstfield2d(fbase,'pbavg   ',pbavg,idm,jdm,0,undef)

      !print *,minval(psikk),maxval(psikk)
      !print *,minval(thkk),maxval(thkk)

      ! cumulative pressure
      p(:,:,:)=0.

      do k=1,kdm

         call read_rstfield2d(fbase,'saln    ',saln ,idm,jdm,k,undef)
         call read_rstfield2d(fbase,'temp    ',temp ,idm,jdm,k,undef)
         call read_rstfield2d(fbase,'dp      ',dp   ,idm,jdm,k,undef)


         ! use upper interface pressure in converting sigma to sigma-star.
         ! this is to avoid density variations in layers intersected by bottom
         ! NB - KAL - thstar/th3d are densities with thbase = 0.
         do j=1,jdm
         do i=1,idm
         if (depths(i,j)>.1 .and. depths(i,j)<1e29) then


            ! sig-option time !
            if (thflag==0) then 
               tmpth3d=   sig0(temp(i,j),saln(i,j))
               tmpkapf=kappaf0(temp(i,j),saln(i,j),p(i,j,k))
            elseif (thflag==2) then 
               tmpth3d=   sig2(temp(i,j),saln(i,j))
               tmpkapf=kappaf2(temp(i,j),saln(i,j),p(i,j,k))
            elseif (thflag==4) then 
               tmpth3d=   sig4(temp(i,j),saln(i,j))
               tmpkapf=kappaf4(temp(i,j),saln(i,j),p(i,j,k))
            else
               print *,'Unknown thflag!!',thflag
               stop '(ssh_from_state)'
            end if

            !th3d(i,j) = sig(temp(i,j),saln(i,j))
            th3d(i,j) = tmpth3d
            if (tbaric) then
               !thstar(i,j,k)=th3d(i,j)+kappaf(temp(i,j), &
               !                       saln(i,j),p(i,j,k))
               thstar(i,j,k)=th3d(i,j)+tmpkapf
            else
               thstar(i,j,k)=th3d(i,j)
            end if
            p(i,j,k+1)=p(i,j,k)+dp(i,j)
            !if (i==200.and.j==300) print *,k,thstar(i,j,k), p(i,j,k+1)/onem
         end if
         enddo
         enddo
      end do

      do j=1,jdm
      do i=1,idm
         if (depths(i,j)>.1 .and. depths(i,j)<1e29) then
         ! store (1+eta) (= p_total/p_prime) in -oneta-
         oneta(i,j)=1.+pbavg(i,j)/p(i,j,kdm+1)

         ! m_prime in lowest layer:
         !montg(i,j,kk)=psikk(i,j)+(p(i,j,kk+1)*(thkk(i,j)-thstar(i,j,kk)) &
         !             -pbavg(i,j)*(thstar(i,j,kk)+thbase))*thref**2
         ! KAL - thkk is density - thbase
         montg(i,j,kdm)=psikk(i,j)+(p(i,j,kdm+1)*(thkk(i,j)+thbase-thstar(i,j,kdm)) &
                      -pbavg(i,j)*(thstar(i,j,kdm)))*thref**2
      end if
      enddo
      enddo

      ! m_prime in remaining layers:
      do k=kdm-1,1,-1
         do j=1,jdm
         do i=1,idm
         if (depths(i,j)>.1 .and. depths(i,j)<1e29) then
            montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*oneta(i,j) &
                    *(thstar(i,j,k+1)-thstar(i,j,k))*thref**2
         end if
         enddo
         enddo
      enddo

      ! SSH
      do  j=1,jdm
      do  i=1,idm
         if (depths(i,j)>.1 .and. depths(i,j)<1e29) then
         ssh(i,j)=(montg(i,j,1)/thref+pbavg(i,j))/onem 
      end if
      end do
      end do

      ! Open unformatted file with ssh
      inquire(iolength=reclSSH) iofld4
      open(20,file='model_SSH.uf',status='unknown', &
           form ='unformatted',access='direct',recl=reclSSH)
      ! Write
      iofld4=ssh
      write(20,rec=imem,iostat=ios2) iofld4

      close(20)
      if(ios2/=0) then
         write(0,'(i1)') 6
         stop '(Could not print to ssh file: p_ssh_from_state.F90)'
      end if
      print *,'field dumped to real*4 binary file model_SSH.uf, record '//cmem
#if defined(DUMP_MSSH)
      open(44,file='lmeanssh.uf',status='replace',form='unformatted')
      iofld=ssh
      write(44) iofld
      close(44)
      print *,'Dump the meanssh in lmeanssh.uf'
#endif
      ! Try to get mean ssh
      if (idm>999 .or. jdm>999) then
         write(tag11,'(i5.5,"x",i5.5)') idm,jdm
      else
         write(tag11,'(i3.3,"x",i3.3)') idm,jdm
      end if
      inquire(exist=ex,file='meanssh'//trim(tag11)//'.uf')
      if(ex) then
         print *,'Found meanssh -- calculating sla'
         open(10,file='meanssh'//trim(tag11)//'.uf',status='old',form='unformatted')
         read(10)iofld
         close(10)
         meanssh=iofld
         print *,minval(iofld),maxval(iofld)
         sla=ssh-meanssh

         ! Open unformatted file with sla
         inquire(iolength=reclSSH)iofld4
         open(20,file='model_SLA.uf',status='unknown', &
              form ='unformatted',access='direct',recl=reclSSH)
         ! Write
         iofld4=sla
         write(20,rec=imem,iostat=ios2) iofld4
         close(20)

         if(ios2/=0) then
            write(0,'(i1)') 6
            stop '(Could not print to sla file: p_ssh_from_state.F90)'
         end if
         print *,'field dumped to real*4 binary file model_SLA.uf, record '//cmem

      else
         print *,'Did not find meanssh --  sla set to zero'
         sla=0.
         meanssh=0.
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TECPLOT STUFF
      open(11,file='tec_SSH.tec',status='unknown')
      write(11,*)'TITLE = "','SSH','"'
      write(11,*)'VARIABLES = "xpos" "ypos" "lon" "lat" "depths" "SSH" "SLA"'
      write(11,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idm,', J=',jdm,', K=1'
 
      write(11,'(30I4)')((          i,i=1,idm),j=1,jdm)
      write(11,'(30I4)')((          j,i=1,idm),j=1,jdm)
      write(11,900)     ((modlon(i,j),i=1,idm),j=1,jdm)
      write(11,900)     ((modlat(i,j),i=1,idm),j=1,jdm)
      write(11,900)     ((depths(i,j),i=1,idm),j=1,jdm)
      write(11,900)     ((ssh   (i,j),i=1,idm),j=1,jdm)
      write(11,900)     ((sla   (i,j),i=1,idm),j=1,jdm)
900   format(11(1x,e12.5))
      close(11)

      print *,'field dumped to tecplot file tec_SSH.tec'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matlab Stuff
!
#if defined (MATLAB)
      ! matlab file
      if (imem==1) then
         matfile='ssh.mat'
      else
         matfile='ssh'//cmem//'.mat'
      end if
      mp=matopen(trim(matfile),'w')
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=modlon;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'modlon', pa1)
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=modlat;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'modlat', pa1)
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=depths;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'depths', pa1)
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=ssh;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'ssh', pa1)
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=meanssh;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'meanssh', pa1)
      !
      pa1=mxCreateNumericMatrix(idm,jdm,mxClassIDFromClassName('double'),0)
      iofld=sla;
      call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),idm*jdm)
      iret=matPutVariable(mp, 'sla', pa1)
      !
      iret=matclose(mp)

      print *,'field dumped to matlab  file '//trim(matfile)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Netcdf Stuff
!
      ! Netcdf - safest bet
      where (depths<.1 .or. depths > 1e20) 
         ssh=undef
         meanssh=undef
         sla=undef
      end where

      ! Undef depth matrix where tiles have been
      where ( depths > 1e20) 
         depths=undef
      end where
      ncfile='sshfromstate.nc'
      print *,'Dumping to netcdf file ',trim(ncfile)
      if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if
      ierr=NF90_DEF_DIM(ncid,'idm',idm,idmid)
      ierr=NF90_DEF_DIM(ncid,'jdm',jdm,jdmid)

      iofld4=modlon
      ierr=NF90_DEF_VAR(ncid,'lon',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4)
 
      iofld4=modlat
      ierr=NF90_REDEF(ncid) ; 
      ierr=NF90_DEF_VAR(ncid,'lat',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 
 
      iofld4=depths
      ierr=NF90_REDEF(ncid) ; 
      ierr=NF90_DEF_VAR(ncid,'depths',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4)) ; 
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=ssh
      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'ssh',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,ssh)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'meanssh',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,meanssh)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'sla',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,sla)

      ierr=NF90_CLOSE(ncid)



      end program p_ssh_from_state
