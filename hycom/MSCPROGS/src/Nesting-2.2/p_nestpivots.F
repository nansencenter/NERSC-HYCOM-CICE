C --- Calculates where the outer model will interpolate its 
c --- model fields onto the grid of the inner model. Inner model is
c --- specified through the nestname variable
      program nestpivots
      use mod_xc 
      use mod_za
      use mod_nest_io
      implicit none
c
      character(len=200) nestname
      integer nrgrids, irec
      logical ex,ass
      integer i,j,k,m
      logical, dimension(:,:), allocatable :: gmsk
      integer, dimension(:,:), allocatable :: imsk
      real   , dimension(:,:), allocatable :: 
     &   plon, plat, depths
      character(len=80) :: cline
      real :: hmin, hmax, hminb, hmaxb
#if defined (IARGC)
      integer*4, external :: iargc
#endif

      if (iargc()==1) then
         call getarg(1,nestname)
      else
         print *,'nestpivots is used to set up the outer model '
         print *,'for nesting. It reads the grid of the outer  '
         print *,'model, then the grid of the inner (nested)   '
         print *,'model. It then calculates the bilinear weights'
         print *,'needed to interpolate the outer model onto   '
         print *,'the inner model grid. This routine should be '
         print *,'run prior to running hycom with outer nesting'
         print *,'enabled. Input argument is path to location  '
         print *,'of inner model regional files. regional files '
         print *,'of outer model must be present in work directory,'
         print *,'as well as the grid.info file for the outer '
         print *,'model'
         print *
         print *,'Usage: nestpivots path_to_inner_model_grid_files'
         call exit(1)
      end if


c --- Read grids of this model
      call xcspmd()
      call zaiost()
      allocate(depths(idm,jdm))
      allocate(plon(idm,jdm))
      allocate(plat(idm,jdm))
      allocate(gmsk(idm,jdm))
      allocate(imsk(idm,jdm))
c
c --- Open and read bathymetry 
      call zaiopf('regional.depth.a', 'old', 701)
      open (unit=701,file='regional.depth.b',status='old')
      call zaiord(depths,imsk,.false.,hmin,hmax,701)
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
          print *, '(forfun_nersc:depth)'
         call exit(1)
       endif
      call zaiocl(701)
      close(701)
c --- Remove 1e30 points
      where(depths>0.5*2.0**99) depths=0.
      gmsk=depths>0.

c --- Open and read p-cell lon and lat
c --- plon
      call zaiopf('regional.grid.a', 'old', 701)
      open (unit=701,file='regional.grid.b',status='old')
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      read(701,'(a80)') cline
      call zaiord(plon,imsk,.false.,hmin,hmax,701)
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
c --- plon
      call zaiord(plat,imsk,.false.,hmin,hmax,701)
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
      call zaiocl(701)
      close(701)
      where (plon<0.) plon=plon+360.


c
c --- TODO: nestname as input argument?
      !nestname='/work/knutali/hycom/FRMa0.05/expt_01.0/SCRATCH/nest/'
c --- Read inner grid also allocates inner_lon/lat/depths
      call readinnergrid(trim(nestname))
c --- Assign pivot points, also allocates:
c --- inner_a1/a2/a3/a4/ipiv/jpiv
      call  assignpivots(trim(nestname), gmsk, plon, plat)

      end program
