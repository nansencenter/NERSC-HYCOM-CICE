#if defined MP_TEST_DUMP
#define MP_TEST_DUMP_CALL call mptestdump
#else
#define MP_TEST_DUMP_CALL !!! call mptestdump(
#endif

!=======================================================================
!BOP
!
! !IROUTINE: evp - elastic-viscous-plastic dynamics driver
!
! !INTERFACE:
!
      subroutine evp(kstrngth)
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic dynamics driver
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use mod_evp
#if defined MP_TEST_DUMP
      use m_mptestdump
#endif
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent(in) ::
     &   kstrngth    
!
!EOP
!
      integer :: k
      integer :: j
      character(len=3) :: ckk

Cdiag ! Test vars
Cdiag real :: hmin, hmax
Cdiag real, dimension(itdm,jtdm) :: guvel, gvvel, guvel2,gvvel2
Cdiag integer, save :: numcalls=0

      call evp_prep(kstrngth)  ! preparation for dynamics

      ! KAL - at this stage the halo should be filled up to nbdy

#if defined MP_TEST_DUMP
      MP_TEST_DUMP_CALL ('evp000' , 0 , 1)
#endif
      imargin=0
      do k=1,ndte              ! subcycling


        write(ckk,'(i3.3)') k

        if (imargin <= 4 ) then
           imargin=nbdy-1 ! Keep open for difference "i-1" in "stress"
           call xctilr(uvel ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
           call xctilr(vvel ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
        end if

#if defined(MIZ)
        call evp_stress_miz(k)     ! the stress tensor equation including MIZ
#elif defined (MIZ2)
        call evp_stress_miz_2(k)   ! the stress tensor equation to model MIZ,
                                   ! but with just shear viscosity
                                   ! reduced in EVP
#else
        call evp_stress(k)         ! the normal EVP stress tensor equation
#endif

        call evp_stepu             ! total surface stress, momentum equation

        imargin=imargin-1

C       MP_TEST_DUMP_CALL ('evp'//ckk , k , 1)


      enddo
#if defined MP_TEST_DUMP
      MP_TEST_DUMP_CALL ('evpndt' , 0 , 1)
#endif

C ============================================================================
C MPI tests
Cdiag numcalls=numcalls+1
Cdiag   if (ijqr == 1) then
Cdiag      call xcaget(guvel,uvel,0)
Cdiag      call xcaget(gvvel,vvel,0)
Cdiag      if (mnproc==1) then 
Cdiag         inquire(iolength=j) guvel,gvvel
Cdiag         open(10,file='evp_tst.uf',status='unknown',
Cdiag&                access='direct',recl=j)
Cdiag         write(10,rec=numcalls) guvel, gvvel
Cdiag         close(10)
Cdiag      end if
Cdiag   else
Cdiag      call xcaget(guvel2,uvel,0)
Cdiag      call xcaget(gvvel2,vvel,0)
Cdiag      if (mnproc==1) then 
Cdiag         inquire(iolength=j) guvel,gvvel
Cdiag         open(10,file='evp_tst.uf',status='old',
Cdiag&             access='direct',recl=j)
Cdiag         read(10,rec=numcalls) guvel, gvvel
Cdiag         close(10)
Cdiag         print *,'****************************'
Cdiag         print *,'EVP diff '
Cdiag         print *,maxval(abs(guvel-guvel2))
Cdiag         print *,maxval(abs(guvel-guvel2))
Cdiag         print *,'****************************'
Cdiag      end if
Cdiag   end if
C ============================================================================
Cdiag call xcstop('evp:mpi-test')

      ! Final tile after everything is finished
      call xctilr(uvel ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(vvel ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      imargin=nbdy-1

      call evp_finish          ! ice-ocean stress

Cdiag print *,'max min uvel ',maxval(uvel ),minval(uvel ),mnproc
Cdiag print *,'max min vvel ',maxval(vvel ),minval(vvel ),mnproc


      

      end subroutine evp
