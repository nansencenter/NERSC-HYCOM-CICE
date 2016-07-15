      module m_mptestdump
      ! This routine contains many variables which use a lot of 
      ! memory. Ex: Topaz (400x600x22) uses 40 MB for one 3D field below
      ! This module is (almost) empty when MP_TEST_DUMP is undefined,
      ! which is usually the case. It only retains the logical parameter
      ! below, which can be used for checking wether the main code is
      ! actually compiled in...
      character(len=30) fnamediff
      character(len=18) fname

      ! KAL -- 20070215 -- Modified to work properly on real*8 archs as well.
      contains
      subroutine mptestdump(tag6,n,MP_TEST_STEP)
      use mod_xc
      use mod_za
      use mod_common_ice
#if defined(ICESTATE)
      use mod_icestate , only : icestate, nthick
#endif
      use mod_forcing_nersc
#if defined(EVP)
      use mod_evp, only: evp_aice=>aice, vice, uocn, vocn, 
     &  uair, vair, ss_tltx, ss_tlty, uvel, vvel,
     &  umassdtei, aiu, strairx, strairy, forcex, forcey,
     &  strintx, strinty,  uarear,
     &  stressp_1, stressp_2, stressp_3, stressp_4,
     &  stressm_1, stressm_2, stressm_3, stressm_4,
     &  stress12_1, stress12_2, stress12_3, stress12_4,
     &  icetmask, iceumask, dxt, dyt, cxm, cym ,cxp, cyp,
     &  HTE, HTN
#endif
#if defined (NOR05) 
      use mod_necessary_ecovars
#endif
      implicit none
      character(len=6),intent(in) :: tag6
      integer, intent(in) :: n,MP_TEST_STEP
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp, tmpdp, 
     &   tmp2
      character(len=5) tag5
      character(len=2) ckk,cll
      character(len=19) fnameout
      logical ex
#if ! defined (MPI) && ! defined (SHMEM)
      integer, external :: omp_get_max_threads
#endif
      integer :: k,l,nop,num_proc,nopdiff
      real    :: tol,xmin,xmax,tmpmax

      logical :: lreal4=.false. ! Routines do not return real 4 versions of numbers
      logical :: ldiff=.false.

      include 'common_blocks.h'
#if defined (NOR05) 
!AS      include 'biocom.h'
#endif


      !  Tolerance value for thin layers
      tol =onem

      if (mod(nstep-nstep1,MP_TEST_STEP)/=0) return

      write(tag5,'(i5.5)')nstep-nstep1

      fname(1:3)='ser'
      fname(4:8)=tag5(1:5)
      fname(9:9)='_'
      fname(10:15)=tag6
      fname(16:18)='.uf'

      fnameout(1:15)=fname(1:15)
      fnameout(16:19)='.dat'

      fnamediff=''
      fnamediff(1:15)=fname(1:15)
      fnamediff=fnamediff//'diff'


      nop = 666
      nopdiff = 667

#if ! defined (MPI) && ! defined (SHMEM)
      num_proc = omp_get_max_threads() 
#else
      num_proc = ijqr
#endif

      inquire(file=fname,exist=ex)
      ! Create file if running on one thread, on one tile  and file doesnt
      ! exist
      !if ((.not.ex).and.ijqr==1) then

      if (mnproc==1) print *,'mptestdump:ijqr,step,id is ',
     &   ijqr,nstep-nstep1,tag6
      !if (ijpr==1) then
      if (.not.ex) then

         !if (mnproc==1) print *,'start write'

         call zaiopf(trim(fname),'replace',nop)
         call zaiowr(depths,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(depthu,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(depthv,ip,.false.,xmin,xmax, nop, lreal4)


         ! Gather data from tiles
         do k=1,4
           call zaiowr(taux(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(tauy(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(uwind(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(vwind(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(airtmp(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(rivers(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
         end do
         call zaiowr(surflx,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(salflx,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(hicem ,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(ficem ,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(hsnwm ,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(iceu  ,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(icev  ,ip,.false.,xmin,xmax, nop, lreal4)
         !if (mnproc==1) print *,' write 1'
#if defined(ICESTATE)
         do k=1,nthick
            tmp=icestate(:,:)%ice(k)%fice
            call zaiowr(tmp,ip,.false.,xmin,xmax, nop, lreal4)
            tmp=icestate(:,:)%ice(k)%hice
            call zaiowr(tmp,ip,.false.,xmin,xmax, nop, lreal4)
            tmp=icestate(:,:)%ice(k)%hsnw
            call zaiowr(tmp,ip,.false.,xmin,xmax, nop, lreal4)
         end do
#endif

#if defined(EVP_MPI)
         call zaiowr(uvel,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(vvel,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(uocn,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(vocn,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(uair,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(vair,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(evp_aice,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(vice,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(ss_tltx,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(ss_tlty,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(umassdtei,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(aiu,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(strairx,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(strairy,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(strintx,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(strinty,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(forcex,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(forcey,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(uarear,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressp_1,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressp_2,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressp_3,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressp_4,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(stressm_1,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressm_2,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressm_3,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stressm_4,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(stress12_1,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stress12_2,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stress12_3,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(stress12_4,ip,.false.,xmin,xmax, nop, lreal4)

         tmp=0. ; where(icetmask) tmp=1.
         call zaiowr(tmp,ip,.false.,xmin,xmax, nop, lreal4)
         tmp=0. ; where(iceumask) tmp=1.
         call zaiowr(tmp,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(dxt,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(dyt,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(cxm,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(cym,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(cxp,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(cyp,ip,.false.,xmin,xmax, nop, lreal4)

         call zaiowr(HTE,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(HTN,ip,.false.,xmin,xmax, nop, lreal4)


         call zaiowr(ULAT,ip,.false.,xmin,xmax, nop, lreal4)
#endif
         call zaiowr(tauxice,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(tauyice,ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(pbavg(:,:,n),ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(ubavg(:,:,n),ip,.false.,xmin,xmax, nop, lreal4)
         call zaiowr(vbavg(:,:,n),ip,.false.,xmin,xmax, nop, lreal4)
         do k=1,kdm
           call zaiowr(dp(:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr( p(:,:,k+1),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(dpu(:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(dpv(:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(u (:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(v (:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(uflx(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(vflx(:,:,k),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(temp(:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
           call zaiowr(saln(:,:,k,n),ip,.false.,xmin,xmax, nop, lreal4)
         end do
#if defined (NOR05)
         if (mnproc==1) print *,'mp_test_dump: dumping ecosys vars'
          do k=1,kdm
          do l=1,nbio
           call zaiowr(bio(:,:,k,n,l),ip,.false.,xmin,xmax, nop, lreal4)
          end do
          end do


          ! KAL - added these "things" as well
!AS          do k=1,kdm
!AS           call zaiowr(prod_dia(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS           call zaiowr(prod_fla(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS           call zaiowr( red_dia(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS           call zaiowr( red_fla(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS           call zaiowr(     rad(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS          end do
!AS          call zaiowr(death_dia(:,:),ip,.false.,xmin,xmax,nop,lreal4)
!AS          call zaiowr(death_fla(:,:),ip,.false.,xmin,xmax,nop,lreal4)
!AS          call zaiowr(diasurf(:,:),ip,.false.,xmin,xmax,nop,lreal4)
!AS          call zaiowr(flasurf(:,:),ip,.false.,xmin,xmax,nop,lreal4)
!AS          do k=1,kdm
!AS             call zaiowr(uflxb(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS             call zaiowr(vflxb(:,:,k),ip,.false.,xmin,xmax,nop,lreal4)
!AS          end do
!AS          call zaiowr(ufluxb,ip,.false.,xmin,xmax,nop,lreal4)
!AS          call zaiowr(vfluxb,ip,.false.,xmin,xmax,nop,lreal4)
#endif /* NOR05 */
             



         !if (mnproc==1) print *,' write 2'

         call zaiocl(nop)
         !if (mnproc==1) print *,' write 3'

      elseif (ex.and.ijpr>1) then
         
         ! Flag if there are differences
         ldiff=.false.


         if (mnproc==1) open(nop,file=fnameout) ! Difference file
         call zaiopf(trim(fname),'old',nop)     ! input file

         if (mnproc==1) then
            write(nop,*) 'Assorted Ocean 2D vars'
            write(nop,*)
         end if

        call compare_flds(depths,'depths',nop,nopdiff,ldiff)
        !call zaiord(tmp,ip,.false.,xmin,xmax, nop)
        !tmp =abs(tmp-real(depths,kind=4))
        !tmpmax=maxval((tmp(1:ii,1:jj)))
        !call xcmaxr(tmpmax)
        !if (mnproc==1) write(nop,*) 'diff depths= ',tmpmax

        call zaiord(tmp,ip,.false.,xmin,xmax, nop)
        tmp =abs(tmp-real(depthu,kind=4))
        tmpmax=maxval((tmp(1:ii,1:jj)))
        call xcmaxr(tmpmax)
        if (mnproc==1) write(nop,*) 'diff depthu= ',tmpmax

        call zaiord(tmp,ip,.false.,xmin,xmax, nop)
        tmp =abs(tmp-real(depthv,kind=4))
        tmpmax=maxval((tmp(1:ii,1:jj)))
        call xcmaxr(tmpmax)
        if (mnproc==1) write(nop,*) 'diff depthv= ',tmpmax

         do k=1,4

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(taux(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff taux= ',tmpmax
           if (abs(tmpmax)/=0.) then
              call zaiowr(tmp,ip,.false.,xmin,xmax, nopdiff, lreal4)
           end if


           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(tauy(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff tauy= ',tmpmax

           call compare_flds(uwind(:,:,k),'uwind',nop,nopdiff,ldiff)
           call compare_flds(vwind(:,:,k),'vwind',nop,nopdiff,ldiff)

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(airtmp(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff airtmp= ',tmpmax

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(rivers(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff rivers= ',tmpmax

         end do

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(surflx(:,:),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff surflx= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(salflx(:,:),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff salflx= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(hicem(:,:),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff hicem= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(ficem(:,:),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff ficem= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(hsnwm(:,:),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff hsnwm= ',tmpmax

         call compare_flds(iceu,'iceu',nop,nopdiff,ldiff)
         call compare_flds(icev,'icev',nop,nopdiff,ldiff)

#if defined(ICESTATE)
         if (mnproc==1) then
            write(nop,*) 'ICESTATE vars'
         end if
         do k=1,nthick

            if (mnproc==1) then
               write(nop,*) 'k=',k
            end if

            call zaiord(tmp,ip,.false.,xmin,xmax, nop)
            tmp2=icestate(:,:)%ice(k)%fice
            tmp=abs(tmp-real(tmp2(:,:),kind=4))
            tmpmax=maxval((tmp(1:ii,1:jj)))
            call xcmaxr(tmpmax)
            if (mnproc==1) write(nop,*) 'diff icestate fice= ',tmpmax

            call zaiord(tmp,ip,.false.,xmin,xmax, nop)
            tmp2=icestate(:,:)%ice(k)%hice
            tmp=abs(tmp-real(tmp2(:,:),kind=4))
            tmpmax=maxval((tmp(1:ii,1:jj)))
            call xcmaxr(tmpmax)
            if (mnproc==1) write(nop,*) 'diff icestate hice= ',tmpmax

            call zaiord(tmp,ip,.false.,xmin,xmax, nop)
            tmp2=icestate(:,:)%ice(k)%hsnw
            tmp =abs(tmp-real(tmp2(:,:),kind=4))
            tmpmax=maxval((tmp(1:ii,1:jj)))
            call xcmaxr(tmpmax)
            if (mnproc==1) write(nop,*) 'diff icestate hsnw= ',tmpmax


            if (mnproc==1) write (nop,*)
         end do
#endif
#if defined(EVP_MPI)
         if (mnproc==1) then
            write(nop,*)
            write(nop,*) 'EVP_MPI vars:'
         end if
         call compare_flds(uvel,'evp_uvel',nop,nopdiff,ldiff)
         call compare_flds(vvel,'evp_vvel',nop,nopdiff,ldiff)
         call compare_flds(uocn,'evp_uocn',nop,nopdiff,ldiff)
         call compare_flds(vocn,'evp_vocn',nop,nopdiff,ldiff)
         call compare_flds(uair,'evp_uair',nop,nopdiff,ldiff)
         call compare_flds(vair,'evp_vair',nop,nopdiff,ldiff)
         call compare_flds(evp_aice,'evp_aice',nop,nopdiff,ldiff)
         call compare_flds(vice,'evp_vice',nop,nopdiff,ldiff)
         call compare_flds(ss_tltx,'evp_sstx',nop,nopdiff,ldiff)
         call compare_flds(ss_tlty,'evp_ssty',nop,nopdiff,ldiff)
         call compare_flds(umassdtei,'evp_umdti',nop,nopdiff,ldiff)
         call compare_flds(aiu,'evp_aiu',nop,nopdiff,ldiff)
         call compare_flds(strairx,'evp_strairx',nop,nopdiff,ldiff)
         call compare_flds(strairy,'evp_strairy',nop,nopdiff,ldiff)
         call compare_flds(strintx,'evp_strintx',nop,nopdiff,ldiff)
         call compare_flds(strinty,'evp_strinty',nop,nopdiff,ldiff)
         call compare_flds(forcex,'evp_forcex',nop,nopdiff,ldiff)
         call compare_flds(forcey,'evp_forcey',nop,nopdiff,ldiff)
         call compare_flds(uarear,'evp_uarear',nop,nopdiff,ldiff)

         call compare_flds(stressp_1,'evp_strp1',nop,nopdiff,ldiff)
         call compare_flds(stressp_2,'evp_strp2',nop,nopdiff,ldiff)
         call compare_flds(stressp_3,'evp_strp3',nop,nopdiff,ldiff)
         call compare_flds(stressp_4,'evp_strp4',nop,nopdiff,ldiff)

         call compare_flds(stressm_1,'evp_strm1',nop,nopdiff,ldiff)
         call compare_flds(stressm_2,'evp_strm2',nop,nopdiff,ldiff)
         call compare_flds(stressm_3,'evp_strm3',nop,nopdiff,ldiff)
         call compare_flds(stressm_4,'evp_strm4',nop,nopdiff,ldiff)

         call compare_flds(stress12_1,'evp_str121',nop,nopdiff,ldiff)
         call compare_flds(stress12_2,'evp_str122',nop,nopdiff,ldiff)
         call compare_flds(stress12_3,'evp_str123',nop,nopdiff,ldiff)
         call compare_flds(stress12_4,'evp_str124',nop,nopdiff,ldiff)

         tmp=0. ; where(icetmask) tmp=1.
         call compare_flds(tmp,'evp_tmsk  ',nop,nopdiff,ldiff)
         tmp=0. ; where(iceumask) tmp=1.
         call compare_flds(tmp,'evp_umsk  ',nop,nopdiff,ldiff)

         call compare_flds(dxt,'evp_dxt  ',nop,nopdiff,ldiff)
         call compare_flds(dyt,'evp_dyt  ',nop,nopdiff,ldiff)

         call compare_flds(cxm,'evp_cxm  ',nop,nopdiff,ldiff)
         call compare_flds(cym,'evp_cym  ',nop,nopdiff,ldiff)
         call compare_flds(cxp,'evp_cxp  ',nop,nopdiff,ldiff)
         call compare_flds(cyp,'evp_cyp  ',nop,nopdiff,ldiff)

         call compare_flds(HTE,'evp_HTE  ',nop,nopdiff,ldiff)
         call compare_flds(HTN,'evp_HTN  ',nop,nopdiff,ldiff)

         call compare_flds(ULAT,'evp_ULAT  ',nop,nopdiff,ldiff)
#endif

         if (mnproc==1) then
            write(nop,*) 'Assorted Ocean 2D vars'
         end if

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(tauxice,kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff tauxice= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(tauyice,kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff tauyice= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(pbavg(:,:,n),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff pbavg= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(ubavg(:,:,n),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff ubavg= ',tmpmax

         call zaiord(tmp,ip,.false.,xmin,xmax, nop)
         tmp =abs(tmp-real(vbavg(:,:,n),kind=4))
         tmpmax=maxval((tmp(1:ii,1:jj)))
         call xcmaxr(tmpmax)
         if (mnproc==1) write(nop,*) 'diff vbavg= ',tmpmax

         if (mnproc==1) then
            write(nop,*)
            write(nop,*) 'Ocean 3D vars:'
         end if

         do k=1,kdm


           call zaiord(tmpdp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmpdp-real(dp(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff dp= ',tmpmax,k

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(p(:,:,k+1),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff  p= ',tmpmax,k


           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(dpu(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff dpu= ',tmpmax,k

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(dpv(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff dpv= ',tmpmax,k


           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(u(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff u= ',tmpmax,k

C          tmp =abs(tmp - u(:,:,k,n))
C          tmpmax=maxval(tmp(1:ii,1:jj),
C    &        tmpdp(1:ii,1:jj)>tol .and. dp(1:ii,1:jj,k,n)>tol)
C          if (mnproc==1) write(nop,*) 'diff u masked= ',tmpmax

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(v(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff v= ',tmpmax,k

C          tmp =abs(tmp - v(:,:,k,n))
C          tmpmax=maxval(tmp(1:ii,1:jj),
C    &        tmpdp(1:ii,1:jj)>tol .and. dp(1:ii,1:jj,k,n)>tol)
C          if (mnproc==1) write(nop,*) 'diff v masked= ',tmpmax


           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(uflx(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff uflx= ',tmpmax,k

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(vflx(:,:,k),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff vflx= ',tmpmax,k

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(temp(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff temp= ',tmpmax,k

C          tmp =abs(tmp - temp(:,:,k,n))
C          tmpmax=maxval(tmp(1:ii,1:jj),
C    &        tmpdp(1:ii,1:jj)>tol .and. dp(1:ii,1:jj,k,n)>tol)
C          if (mnproc==1) write(nop,*) 'diff temp masked= ',tmpmax

           call zaiord(tmp,ip,.false.,xmin,xmax, nop)
           tmp =abs(tmp-real(saln(:,:,k,n),kind=4))
           tmpmax=maxval((tmp(1:ii,1:jj)))
           call xcmaxr(tmpmax)
           if (mnproc==1) write(nop,*) 'diff saln= ',tmpmax,k

C          tmp =abs(tmp - saln(:,:,k,n))
C          tmpmax=maxval(tmp(1:ii,1:jj),
C    &        tmpdp(1:ii,1:jj)>tol .and. dp(1:ii,1:jj,k,n)>tol)
C          if (mnproc==1) write(nop,*) 'diff saln masked= ',tmpmax

           if (mnproc==1)  write(nop,*)
         end do
#if defined (NOR05)
          do k=1,kdm
          do l=1,nbio
             write(ckk,'(i2.2)') k
             write(cll,'(i2.2)') l
             call compare_flds(bio(:,:,k,n,l),
     &           'bio '//ckk//' '//cll//' :', nop,nopdiff,ldiff)
          end do
          end do
#endif /* NOR05 */

        if (ldiff) then
           call zaiocl(nopdiff)
           if (mnproc==1) close(nopdiff)
        end if
        call zaiocl(nop)
        if (mnproc==1) close(nop)

      endif
      !if (mnproc==1) print *,'end mptestdump'
      end subroutine mptestdump




      subroutine compare_flds(fld,cfld,nop,nopdiff,ldiff)
      use mod_xc
      use mod_za
      implicit none
      character(len=*), intent(in)   :: cfld
      real,             intent(in)   ::
     &   fld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer,          intent(in)   :: nop,nopdiff
      logical,          intent(inout):: ldiff

      real    :: tmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    :: tmpmax,xmin,xmax

      call zaiord(tmp,ip,.false.,xmin,xmax, nop)
      tmp =abs(tmp-real(fld,kind=4))
      ! We are only interested in the ocean. Forcing fields may
      ! cause problems here
      where(ip==0) tmp=0.
      tmpmax=maxval((tmp(1:ii,1:jj)))
      call xcmaxr(tmpmax)
      if (mnproc==1) write(nop,*) 'diff '//trim(cfld)//'= ',tmpmax
      if (abs(tmpmax)/=0.) then


         ! Open in we discovered a difference
         if (.not.ldiff) then
            call zaiopf(trim(fnamediff)//'.a','replace',nopdiff)
            if (mnproc==1)
     &         open(nopdiff,file=trim(fnamediff)//'.b',status='replace')
            ldiff =.true.
         end if

         call zaiowr(tmp,ip,.false.,xmin,xmax, nopdiff, .false.)
         if (mnproc==1) write(nopdiff,'(a)') cfld
      end if
      end subroutine compare_flds


      end module m_mptestdump
