      program nestini
C --- Nestini sets up the relaxation factors 
C --- which are used in nestbnd and nestbarotp
C ---
C --- This routine works on the global grid -- always
      use mod_xc
      use mod_za
      use netcdf
      use m_parse_blkdat
      implicit none
      integer i,j,i2,j2,l
      real, allocatable, dimension(:,:) :: depths,
     &   nestrel, nestrelu, nestrelv, factor
      integer, allocatable, dimension(:,:) ::  idone
      logical, allocatable, dimension(:,:) ::  ldone
      real , allocatable :: nestfunc(:)
      integer :: ios,iyear,iday, inest
      logical :: ex
      character(len=80) :: cfile
      integer :: ncid, varid, idmid,jdmid
      real    :: hmin, hmax, hminb, hmaxb
      character(len=80) :: cline
      real, parameter :: undef=-1e14
      real    :: relax_strength=0.90 !Fanf: add the possibillity to tune
      !the relax strength with a non dimensional parameter (1 means no impact)
      integer :: idummy
      real    :: nestfq=0.25 !in day
c
c --- Init io and arrays
      inest=20
      call xcspmd()
      call zaiost()
      allocate(depths(idm,jdm))
      allocate(nestfunc(inest))
      allocate(nestrel (idm,jdm))
      allocate(nestrelu(idm,jdm))
      allocate(nestrelv(idm,jdm))
      allocate(factor  (idm,jdm))
      allocate(idone(idm,jdm))
      allocate(ldone(idm,jdm))
      nestrel=0.0
c
c --- Open and read bathymetry
      call zaiopf('regional.depth.a', 'old', 701)
      open (unit=701,file='regional.depth.b',status='old')
      call zaiord(depths,idone,.false.,hmin,hmax,701)
      read(701,'(a80)') cline ; print *,cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      i = index(cline,'=')
      read (cline(i+1:),*) hminb,hmaxb
      if (abs(hmin-hminb).gt.abs(hminb)*1.e-4 .or. 
     &    abs(hmax-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') 
     &      'error - .a and .b files not consistent:', 
     &      '.a,.b min = ',hmin,hminb,hmin-hminb, 
     &      '.a,.b max = ',hmax,hmaxb,hmax-hmaxb
          print *, '(forfun_nersc:plon)'
         call exit(1)
       endif
      ! Remove 1e30 points
      call zaiocl(701)
      close(701)
      where (depths>1e20) depths=0.

c      Get the baclin
      !call parse_blkdat('baclin','integer',rvar,baclin)
       call parse_blkdat('nestfq','real',nestfq,idummy)
       if (nestfq .lt. 0.01) then ! Unlikely high freq 
         print *,'nestfq in blkdat.input is 0, quit'
        call exit(1)
       endif
      
c
c --- relaxation factor -- depends on distance from boundary
      do i=1,inest
         nestfunc(i)=1.-float(i-1)/float(inest-1)
         !nestfunc(i)=min(1.0,max(0.,nestfunc(i))**2)
         ! Yves
         nestfunc(i)=relax_strength*min(1.0,max(0.,nestfunc(i))**4)

         !Fanf: We multiply here to ensure the follwing formulae
         ! P_inner=P_inner+trelax*baclin/(nestfq in sec)(P_inner-P_outer)
         !note that in termf we have :
         ! P_inner=P_inner+2*baclini*nestfunc(P_inner-P_outer)
         nestfunc(i)=nestfunc(i)*1./(2.*nestfq*86400)
      enddo
c
c
c --- Relaxation functions given values dependent on the 
c --- distance to the closes open boundary point
      do j=2,jdm-1
         if (depths(2,j) > 1.0) then 
            nestrel(2,j)=nestfunc(1)
            idone(2,j)=1
            ldone(2,j)=.true.
         endif
         if (depths(idm-1,j) > 1.0) then
            nestrel(idm-1,j)=nestfunc(1)
            idone(idm-1,j)=1
            ldone(idm-1,j)=.true.
         endif
      enddo
      do i=2,idm
         if (depths(i,2) > 1.0) then
            nestrel(i,2)=nestfunc(1)
            idone(i,2)=1
            ldone(i,2)=.true.
         endif
         if (depths(i,jdm-1) > 1.0) then
            nestrel(i,jdm-1)=nestfunc(1)
            idone(i,jdm-1)=1
            ldone(i,jdm-1)=.true.
         endif
      enddo
c
c --- nest factor in p-points
      do l=2,inest-1
         do j=2,jdm-1
         do i=2,idm-1
            if (idone(i,j)==l-1) then
              i2=max(i-1,2)
              j2=j   
              if (depths(i2,j2) > 0.0 .and. .not.ldone(i2,j2)) then
                 ldone(i2,j2)=.true.
                 idone(i2,j2)=l
                 nestrel(i2,j2)=nestfunc(l)
              endif

              i2=min(i+1,idm-1)
              j2=j   
              if (depths(i2,j2) > 0.0 .and.  .not.ldone(i2,j2)) then
                 ldone(i2,j2)=.true.
                 idone(i2,j2)=l
                 nestrel(i2,j2)=nestfunc(l)
              endif

              i2=i
              j2=max(j-1,2)
              if (depths(i2,j2) > 0.0  .and. .not.ldone(i2,j2)) then
                 ldone(i2,j2)=.true.
                 idone(i2,j2)=l
                 nestrel(i2,j2)=nestfunc(l)
              endif

              i2=i
              j2=min(j+1,jdm-1)
              if (depths(i2,j2) > 0.0 .and.  .not.ldone(i2,j2)) then
                 ldone(i2,j2)=.true.
                 idone(i2,j2)=l
                 nestrel(i2,j2)=nestfunc(l)
              endif

            endif
         enddo
         enddo
      enddo
c
c --- Set up nesting factor in u points
      do i=2,idm
         nestrelu(i,:)=0.25*(nestrel(i,:)+nestrel(i-1,:))
      enddo
      nestrelu(1,:)=0.0
      nestrelu(:,1)=0.0
      nestrelu(:,jdm)=0.0
c
c --- Set up nesting factor in v points
      do j=2,jdm
         nestrelv(:,j)=0.25*(nestrel(:,j)+nestrel(:,j-1))
      enddo
      nestrelv(1,:)=0.0
      nestrelv(:,1)=0.0
      nestrelv(idm,:)=0.0
c
c --- Set to zero at zero depths
      where (depths(2:idm-1,2:jdm-1) == 0.0)
         nestrelu(2:idm-1,2:jdm-1)=0.0
         nestrelv(2:idm-1,2:jdm-1)=0.0
      endwhere
c
      cline='Single relaxation mask'
      call zaiopf('rmu.a', 'replace', 915)
      open (unit=915,file='rmu.b',status='replace')
      write(915,'(a80)') cline
      write(915,'(a80)') 
      write(915,'(a80)') 
      write(915,'(a80)') 
      write(915,'(a80)') 
      call zaiowr(nestrel,idone,.false.,hmin,hmax,915,.false.)
      write(cline,'("rmu =",2f24.12)') hmin, hmax
      write(915,'(a80)') cline
      close(915)
      call zaiocl(915)
c
      factor=0.
      where (ldone) factor=1.0
      where (depths<.1) nestrel=undef
      where (depths<.1) nestrelu=undef
      where (depths<.1) nestrelv=undef
      where (depths<.1) factor=undef
c
c --- Netcdf Diagnostics
      call ncerr(NF90_create('nestrel.nc',NF90_CLOBBER,ncid))
      call ncerr(NF90_def_dim(ncid,'idm',idm,idmid))
      call ncerr(NF90_def_dim(ncid,'jdm',jdm,jdmid))
      call ncerr(NF90_def_var(ncid,'depths',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',0.))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,depths))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'nestrel',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,nestrel))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'nestrelu',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,nestrelu))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'nestrelv',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,nestrelv))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'factor',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,factor))
      call ncerr(NF90_close(ncid))
c
      print *,"Finished setting up nesting relaxation factors rmu.[ab]"
      print *,"and diagnostics file nestrel.nc"

c
      contains
         subroutine ncerr(errcode)
         use netcdf
         implicit none
         integer, intent(in) :: errcode
         if (errcode/=NF90_NOERR) then
            write(lp,'(a)') NF90_STRERROR(errcode)
            stop '(ncerr)'
         end if
         end subroutine
      end program
