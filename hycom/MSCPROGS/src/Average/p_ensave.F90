program p_ensave
   !------------------------------------------------
   !Routine ensave:
   !---------------
   !Routine creates a file containing mean values 
   !of the variables contained in the input files. 
   !All input files must be of the same type.  Note 
   !that the means of layer (3D ocean) variables 
   !are weighted by layer thickness.
   !
   ! Created by Francois Counillon, 2012
   !------------------------------------------------
   use mod_xc
   use mod_za
   use mod_parameters
   use mod_hycomfile_io
   implicit none
   character(len=80) :: infile,ftype,outfile
#if defined (IARGC)
   integer*4, external :: iargc
#endif
   integer :: nfiles,ifile
   real, dimension(:,:), allocatable :: field,flag,mean,dp, sumdp
   real*8, dimension(:,:), allocatable :: iofld
   character(len=8) :: cfld
   character(len=3) :: cmem,ctmp
   integer          :: coord, ierr, indx, nstep, tlevel, ensize,imem
   real   :: dens, rday, adens
   type(hycomfile) :: df0, df, dfave

   if (iargc()==0) then
      print *,'Routine creates a file containing mean values of the '
      print *,'variables contained in the input files. All input files'
      print *,'must be of the same type. All date must be the same'
      print *,'Note that the means of layer variables are weighted by'
      print *,'Layer thickness'
      print *,'Usage:'
      print *,'ensave filetype file1 ensize'
      print *
      print *,'filetype is either restart, nersc_daily or nersc_weekly'
      print *,'Example:'
      print *,'   ensave nersc_daily TP3DAILY_1990_240_1990_313 100 '
      print *
      call exit(1)
    end if

   call xcspmd()
   call zaiost()

   allocate(field(idm,jdm))
   allocate(dp   (idm,jdm))
   allocate(sumdp   (idm,jdm))
   allocate(flag (idm,jdm))
   allocate(mean (idm,jdm))
   ! init output file to something .. unlikely ..
   call getarg(1,ftype)
   call getarg(2,outfile)
   call getarg(3,ctmp)
   read(ctmp,*) ensize

   if (trim(ftype)/='nersc_weekly' .and. trim(ftype)/='nersc_daily' .and. trim(ftype)/='restart') then
      print *,'ftype'
      print *,'Routine creates a file containing mean values of the '
      print *,' an ensemble of input files. All file must be at the same date'
      print *,'must be of the same type.'
      print *,'Note that the means of layer variables are weighted by'
      print *,'Layer thickness'
      print *,'Usage:'
      print *,'ensave filetype filebase ensize'
      print *
      print *,'filetype is either restart, nersc_daily or nersc_weekly'
      print *,'Example:'
      print *,'   ensave nersc_daily TP3DAILY_1990_240_1990_313 100'
      print *
      call exit(1)
    end if


   ! First file is used for reading fields, and for initializing the averaged file 
   call initHF(df0,trim(outfile)//'_mem001.a',trim(ftype))
   call initHF(dfave,trim(outfile)//'.a',trim(ftype),template=df0) 

   ierr=0
   indx=1
   do while (ierr==0) 

      ! Get field to process - use first file as template
      call HFHeaderFromIndex(df0,indx,ierr,varname=cfld,level=coord,timelevel=tlevel)
      if (ierr/=0) then
         print *,'---'
         print *,'Normal exit - reached end of first input file'
         print *,'Average file is : '//trim(outfile)
         call exit(0)
      end if
      ! process mean of this field for all files
      mean  = 0.
      sumdp = 0.
      ifile=1
      do imem=1,ensize
         write(cmem,'(i3.3)') imem
         call initHF(df,trim(outfile)//'_mem'//cmem//'.a',ftype)
         call HFReadField(df,field,idm,jdm,cfld,coord,tlevel)
         !print *,cfld,minval(field),maxval(field)
         !print '("---->",a8,3i5,2e14.4)',cfld,coord,indx,ierr,minval(field),maxval(field)

         ! Flag undefined values
         flag=1. ;  where (field>1e28) flag=0.

         ! Accumulate layer-weighted/arithmetic sum
         if (is3dvar(df,cfld,tlevel)   .and. .not. isDPVar(df,cfld)  ) then
            call HFReadDPField(df,dp,idm,jdm,coord,tlevel)
            ! More robust mean wrt precision loss
            mean =mean +(field*dp*flag-mean)/ifile
            sumdp=sumdp+(dp-sumdp*flag)/ifile
         !   print *,'toto',mean(13,137),sumdp(13,137),field(13,137),dp(13,137)
         else
            ! More robust mean wrt precision loss
            mean=mean+(field*flag-mean)/ifile
            sumdp=1.
         end if

         ! On first pass, do consistency check on files
         if (indx==1) then
            call HFUpdateAverage(dfave,df)
         end if 
         ifile=ifile+1
      end do ! ens loop
      ! On first pass , dump FILE header
      if (indx==1) then
         call HFWriteHeader(dfave,vDim(df))
      end if

      ! Calculate layer-weighted/arithmetic mean
      if (is3dvar(df,cfld,tlevel) .and. .not.isDPVar(df,cfld)) then

         ! nersc_weekly have the layer sum on disk, others have layer average 
         if (trim(ftype)/='nersc_weekly') then
            mean=mean/max(sumdp,1e-4)
         end if
      end if
      ! .. and write field to disk (also writes FIELD header to disk)
      call HFWriteField(dfave,mean,idm,jdm,cfld,coord,tlevel,indx)
      indx=indx+1
      print '(a8,3i5,2e14.4)',cfld,coord,indx,ierr,minval(mean),maxval(mean)
   end do
end program p_ensave
