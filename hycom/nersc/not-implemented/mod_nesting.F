C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  This file contain various routines needed for generating nesting conditions
C  using a "outer" model and the routines needed for applying the
C  nesting conditions in the local high resolution model.
C
C  The following routines are defined here:
C  1.  module mod_nesting:
C      definitions of relaxation variables and flags
C
C  2.  subroutine nestcondread:
C      read nesting conditions from file.  Used by inner model.
C
C  3.  subroutine nestcondsave:
C      Reads locations of nesting grid points, interpolates and saves
C      the nesting variables.  Used by the outer model.
C
C  The above routines are changed in this version so that they are easier
C  to use (meaning they are understandable to other than the original
C  author ;-)
C
C  The main purpose of this module is to be a drop-in replacement for
C  the hycom routines "rdnest_in" and "rdbaro_in" found in forfun.F
C
C  Note that the routines rely heavily on mod_confmap and must be
C  changed if conformal mepping is not used.
C
C  Also note that optimization is probably possible, especially on the
C  output/input routines. So far this seems not to be a huge issue but you
C  may feel different...
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  Auxillary routines:
C  function   filenesting2  - returns name of nesting file to read/save
C  subroutine writenest     - Writes one variable into nesting files
C  subroutine readnest      - Reads  one variable from nesting file
C  subroutine interpnest    - Takes care of interpolation to nesting points
C  subroutine readinnergrid - Reads inner model nesting points
C  subroutine readpivots    - Assigns pivot points for interpolation from
C                            outer model
C
C  Note that these routines are private to the module
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Routines that USED to be here
C
C  subroutine nestloc_save:
C      Dumped the locations of the nesting gridpoints to file.  Used by
C      the local model initially and file is later read by global model.
C      KAL - Now done as a preprocessing step
C
C  subroutine nestbnd:
C      Applied the nesting conditions in the local model.
C      KAL - Now uses the nesting procedure used in standard hycom (see thermf)
C
C  subroutine nestbarotp:
C      Applied the barotropic nesting conditions in the local model.
C      KAL - Now uses the procedure used in standard hycom (see barotp/latbd)
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module mod_nesting
      use mod_xc
      implicit none

      integer, parameter :: inest=20  ! Width of boundary zone (grid points)
      logical, save :: lnesto         ! Save outer nesting bondary conditions
      integer*8, save :: nestdto        ! dt of outer nesting conditions
      logical, save :: lnesti         ! Read and apply nesting bondary conditions
      integer, save :: nestdti        ! dt of inner nesting conditions


c --- Variables used when saving nesting conditions. Convenient to 
c --- keep them here
      integer, save, private,  allocatable, dimension(:,:) :: 
     &   inner_ipiv, inner_jpiv
      real,    save, private,  allocatable, dimension(:,:) :: 
     &   inner_a1, inner_a2, inner_a3, inner_a4, inner_depths,
     &   inner_lon, inner_lat
      logical, save, private, allocatable, dimension(:,:) :: 
     &   inner_tile_mask
      integer, save, private :: inner_idm, inner_jdm, 
     &   inner_iidm, inner_jjdm, inner_inest

      private :: filenesting2, writenest, readnest, interpnest, 
     &           readinnergrid
      contains

c
C --- Save nesting conditions from this model to inner models,
C --- Location of saves specified in infile read below (nestname)
C --- NB - this routine will dump everything, meaning 3D fields and 
C --- Barotropic fields
C
      subroutine nestcondsave(dtime,n,init,imember)
      use mod_xc 
      use mod_za
      use mod_year_info, only: year_info
      use mod_hycom_nersc
#if defined (ICE) && defined (ICE_NEST)
      use mod_common_ice
#endif
      implicit none
      real*8 , intent(in)    :: dtime
      integer, intent(in)    :: n
      logical, intent(in)    :: init
      integer                :: imember
c
      integer, parameter :: max_nrgrids=100
      character(len=80) nestname(max_nrgrids)
      integer nrgrids, irec
      logical ex,ass
      integer ipiv,jpiv,gipiv,gjpiv
      integer i,j,k,m
      integer :: iostatus
#if defined (NEST_BIO)
#if defined (NOR05)
      integer ktr
      character trcnum*2
#endif /* NOR05 */
#endif /* NEST_BIO */
      real, allocatable, dimension(:,:) :: 
     &   inner_tmp1,   inner_tmp1a,
     &   inner_tmp2,   inner_tmp2a
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   fldsum, fld1, fld2
      character(len=200) :: ncfil
      character(len=10) :: ctt
      integer*8 :: ihour,isec
      real*8  :: dtime2
      real :: xmin, xmax
      include 'common_blocks.h'
      include 'stmt_fns.h'
c
c --- Check if this is a nest output time
c --- This has the benefit of working when 
c --- int(24/nestdto)*nestdto /= 24/nestdto
      ihour=nint(dtime*86400.d0,kind=8)/3600
      isec =nint(dtime*86400.d0,kind=8)-ihour*3600
cdiag if (mnproc==1) print *,'nestcondsave hour, sec:',ihour,isec
      if (mod(ihour,nestdto)==0 .and. isec < baclin) then
         dtime2=(ihour*3600.d0)/86400.d0
         continue
C      else if (init) then
C         dtime2=(ihour/nestdto)*real(nestdto,kind=8)
C         continue
      else 
         return
      end if
c
      if (mnproc==1) write(lp,'(a)') 'saving outer nesting conditions '
c
      ! read dimensions and locations for internal grids   
      inquire(file='nesting.in',exist=ex)
      if (.not.ex) stop 'nestcond_save: nesting.in does not exist'
      open(43,file='nesting.in',action='read',STATUS='OLD',
     &  iostat = iostatus)
         do m=1,max_nrgrids
            read(43,'(a)', iostat = iostatus ) nestname(m)
            if (iostatus /= 0) then
               nrgrids=m-1
               exit
            endif
            if (mnproc==1) print '(a,i3)','nestname read from file= '//
     &         trim(nestname(m))//' for grid no ',m
         enddo
      close(43)
c
c --- Loop through the internal grids specified in nestname files (read above)
      do m=1,nrgrids
      inquire(file=trim(nestname(m))//'regional.grid.a',exist=ex)
      if (.not. ex) then
         if (mnproc==1) then
            write(lp,*) 'Could not find'//
     &         trim(nestname(m))//'regional.grid.a'
            write(lp,*) 'I quit'
         end if
         call xcstop('nestcondsave')
         stop '(nestcondsave)'
      else 
c ---    Read inner grid also allocates inner_lon/lat/depths
         call readinnergrid(nestname(m))
c ---    Read preprocessed pivot points, also allocates:
c ---    inner_a1/a2/a3/a4/ipiv/jpiv/tile_mask
         call  readpivots(nestname(m))
c
c ---    Set up temporary arrays for outer data
         allocate(inner_tmp1  (inner_idm,inner_jdm))
         allocate(inner_tmp1a (inner_idm,inner_jdm))
         allocate(inner_tmp2  (inner_idm,inner_jdm))
         allocate(inner_tmp2a (inner_idm,inner_jdm))
c
c ---    -----------------------------------------------------------------
c ---    output 3d nesting fields and 2D (barotropic/ice) fields
c ---    New and more intuitive approach
c ---    -----------------------------------------------------------------
c ---    temp
         irec=1
         do k=1,kdm
            call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &         temp(1-nbdy,1-nbdy,k,n),dp(1-nbdy,1-nbdy,k,n),k>1)
            call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
            call writenest(nestname(m),dtime2,'temp    ',k,
     &         inner_tmp1,inner_idm,inner_jdm,irec,imember)
            irec=irec+1
         end do
c ---    saln
         do k=1,kdm
            call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &         saln(1-nbdy,1-nbdy,k,n),dp(1-nbdy,1-nbdy,k,n),k>1)
            call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
            call writenest(nestname(m),dtime2,'saln    ',k,
     &         inner_tmp1,inner_idm,inner_jdm,irec,imember)
            irec=irec+1
         end do
c ---    u/v (baroclinic)
         call xctilr(u(1-nbdy,1-nbdy,1,n),1,kk, nbdy,nbdy, halo_uv)
         call xctilr(v(1-nbdy,1-nbdy,1,n),1,kk, nbdy,nbdy, halo_vv)
         do k=1,kdm
            fld1=u(:,:,k,n)
            fld2=v(:,:,k,n)

CKAL        if (k==1) then
CKAL           write(ctt,'(f10.2)') dtime
CKAL           call zaiopf('tstubef'//trim(adjustl(ctt))//'.uf',
CKAL &            'replace',999)
CKAL           call zaiowr(fld1,ip,.false.,xmin,xmax, 999, .true.)
CKAL           call zaiocl(999)
CKAL           call zaiopf('tstvbef'//trim(adjustl(ctt))//'.uf',
CKAL &            'replace',999)
CKAL           call zaiowr(fld2,ip,.false.,xmin,xmax, 999, .true.)
CKAL           call zaiocl(999)
CKAL        end if

            call rotate2(fld1,fld2,plat,plon,'m2l',.false.)
            !call rotate(fld1,fld2,plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
            call xctilr(fld1(1-nbdy,1-nbdy),1,1, nbdy,nbdy, halo_uv)
            call xctilr(fld2(1-nbdy,1-nbdy),1,1, nbdy,nbdy, halo_vv)

CKAL        if (k==1) then
CKAL           write(ctt,'(f10.2)') dtime
CKAL           call zaiopf('tstuaft'//trim(adjustl(ctt))//'.uf',
CKAL &            'replace',999)
CKAL           call zaiowr(fld1,ip,.false.,xmin,xmax, 999, .true.)
CKAL           call zaiocl(999)
CKAL           call zaiopf('tstvaft'//trim(adjustl(ctt))//'.uf',
CKAL &            'replace',999)
CKAL           call zaiowr(fld2,ip,.false.,xmin,xmax, 999, .true.)
CKAL           call zaiocl(999)
CKAL        end if

c
            call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &         fld1(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,k,n),k>1)
            call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
            call writenest(nestname(m),dtime2,'u       ',k,
     &         inner_tmp1,inner_idm,inner_jdm,irec,imember)
            irec=irec+1
c
            call interpnest(inner_tmp2,inner_tmp2a,inner_idm,inner_jdm,
     &         fld2(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,k,n),k>1)
            call xcnest_reduce(inner_tmp2,inner_idm,inner_jdm)
            call writenest(nestname(m),dtime2,'v       ',k,
     &         inner_tmp2,inner_idm,inner_jdm,irec,imember)
            irec=irec+1
         end do
c ---    intf
         fldsum=0.
         do k=1,kdm
            fldsum=fldsum+dp(:,:,k,n)
            call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &         fldsum(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,k,n),.false.)
            call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
            call writenest(nestname(m),dtime2,'intf    ',k,
     &         inner_tmp1,inner_idm,inner_jdm,irec,imember)
            irec=irec+1
         end do
c ---    Barotropic velocity
         call xctilr(ubavg(1-nbdy,1-nbdy,n),1,1, nbdy,nbdy, halo_uv)
         call xctilr(vbavg(1-nbdy,1-nbdy,n),1,1, nbdy,nbdy, halo_vv)
         fld1=ubavg(:,:,n)
         fld2=vbavg(:,:,n)
         call rotate(fld1,fld2,plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
c
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      fld1(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'ubavg   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      fld2(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'vbavg   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c ---    Barotropic pressure
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      pbavg(1-nbdy,1-nbdy,n),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'pbavg   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
#if defined (ICE) && defined (ICE_NEST)
c ---    Ice Concentration
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      ficem(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'ficem   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c ---    Ice Thickness
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      hicem(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'hicem   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c ---    Ice Thickness
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      hsnwm(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'hsnwm   ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c ---    Ice velocity
         call xctilr(iceu,1,1, nbdy,nbdy, halo_uv)
         call xctilr(icev,1,1, nbdy,nbdy, halo_vv)
         fld1=iceu
         fld2=icev
         call rotate(fld1,fld2,plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      fld1(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'uice    ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
c
         call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &      fld2(1-nbdy,1-nbdy),dp(1-nbdy,1-nbdy,1,n),.false.)
         call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
         call writenest(nestname(m),dtime2,'vice    ',0,
     &      inner_tmp1,inner_idm,inner_jdm,irec,imember)
         irec=irec+1
#endif
#if defined (NEST_BIO)
#if defined (NOR05)
c ---    tracer
         do ktr=1,ntracr
           write(trcnum,'(i2.2)') ktr
           do k=1,kdm
             call interpnest(inner_tmp1,inner_tmp1a,inner_idm,inner_jdm,
     &         tracer(1-nbdy,1-nbdy,k,n,ktr),dp(1-nbdy,1-nbdy,k,n),k>1)
             call xcnest_reduce(inner_tmp1,inner_idm,inner_jdm)
             call writenest(nestname(m),dtime2,'trcr'//trcnum//'  ',k,
     &           inner_tmp1,inner_idm,inner_jdm,irec)
             irec=irec+1
           end do
         end do
#endif /* NOR05 */
#endif /* NEST_BIO */
c ---    -----------------------------------------------------------------
c ---    Deallocate temporary arrays for next nest output
         deallocate( inner_depths   )
         deallocate( inner_lon      )
         deallocate( inner_lat      )
         deallocate( inner_tmp1     )
         deallocate( inner_tmp1a    )
         deallocate( inner_tmp2     )
         deallocate( inner_tmp2a    )
         deallocate( inner_a1       )
         deallocate( inner_a2       )
         deallocate( inner_a3       )
         deallocate( inner_a4       )
         deallocate( inner_ipiv     )
         deallocate( inner_jpiv     )
         deallocate( inner_tile_mask)
      endif
      enddo ! Next nest grid

      if (mnproc==1) then 
         write(lp,'(a)') 'Finished dumping nesting conditions ..'
      end if
      end subroutine nestcondsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

      subroutine nestcondread(dtime,lslot,btrpflag,bclnflag)
      use mod_xc
#if defined (ICE_NEST) && defined (ICE)
      use mod_common_ice
#endif
#if defined (NEST_BIO)
#if defined (NOR05)
      use mod_necessary_ecovars
#endif /* NOR05 */
#endif /* NEST_BIO */
      use mod_hycom_nersc
      implicit none
      real*8         , intent(in) :: dtime
      integer        , intent(in) :: lslot
      logical        , intent(in) :: btrpflag
      logical        , intent(in) :: bclnflag
c
      character(len=3) :: ckk
      character(len=10) :: ctt
      character(len=80) :: ncfil
      real, dimension(itdm,jtdm) :: fld,fld2
      real :: hqpbot, hqpbotu, hqpbotv
      integer j,k,i
#if defined (NEST_BIO)
#if defined (NOR05)
      integer :: ktr
      character :: trcnum*2
#endif /* NOR05 */
#endif /* NEST_BIO */
      logical :: errflg
      logical, parameter :: diag=.false.
      logical, save      :: first=.true.
      include 'common_blocks.h'
c
      if (first) then
c ---    Error check on nest dimensions
         if (kknest/=kdm) then
            if (mnproc==1) then
               write(lp,'(a)') 'kknest is not=kdm '
               write(lp,'(a)') 'Fix this in order to use nesting '
               call flush(lp)
            end if
            call xcstop('(nestcondread)')
            stop '(nestcondread)'
         end if
c ---    Error check on rmu masks
         call xcaget(fld ,rmunp,0)
         call xcaget(fld2,rmunv,0)
         errflg=.false.
         do j=1,jtdm
         do i=1,itdm
            k=min(itdm-i,i,jtdm-j,j)
            if ((fld(i,j)>1e-16 .or. fld2(i,j)>1e-16).and. k>inest) then
               if (mnproc==1) 
     &         write(lp,'(a,2i5)') 'Nest boundary error at ',
     &         i,j
               errflg=.true.
            end if
         end do
         end do
         if (errflg) then
            if (mnproc==1) write(lp,'(a)') 'Nestmask error '
            call xcstop ('(mod_nesting:nestcondread)')
            stop '(mod_nesting:nestcondread)'
         end if
         first=.false.
      end if
c
      if (mnproc==1) then
         write(lp,'(a,f10.2)') 'Nestcondread called at ',dtime
      end if
c
c --- Read 3d and baroclinic nesting fields
      if (bclnflag) then 
c ---    temp - mnproc=1 returns field
         do k=1,kdm
           call readnest('nest/',dtime,'temp    ',k, fld,itdm,jtdm,imem)
           call xcaput(fld,tnest(1-nbdy,1-nbdy,k,lslot),1)
         end do
c ---    saln - mnproc=1 returns field
         do k=1,kdm
           call readnest('nest/',dtime,'saln    ',k, fld,itdm,jtdm,imem)
           call xcaput(fld,snest(1-nbdy,1-nbdy,k,lslot),1)
         end do
c ---    u - mnproc=1 returns field
         do k=1,kdm
           call readnest('nest/',dtime,'u       ',k, fld,itdm,jtdm,imem)
           call xcaput(fld,unest(1-nbdy,1-nbdy,k,lslot),1)
         end do
c ---    v - mnproc=1 returns field
         do k=1,kdm
           call readnest('nest/',dtime,'v       ',k, fld,itdm,jtdm,imem)
           call xcaput(fld,vnest(1-nbdy,1-nbdy,k,lslot),1)
         end do
c ---    pressure - mnproc=1 returns field
c ---    NB: intf is lowest interface of a layer, but pnest is top
c ---    interface of layer
c ---    TODO: change this behaviour to be more comnsistent with pnest
         pnest(:,:,1,lslot)=0.
         do k=1,kdm-1
           call readnest('nest/',dtime,'intf    ',k, fld,itdm,jtdm,imem)
           call xcaput(fld,pnest(1-nbdy,1-nbdy,k+1,lslot),1)
         end do
#if defined (ICE_NEST) && defined (ICE)
c ---    hicem    - mnproc=1 returns field
         call readnest('nest/',dtime,'hicem   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,hicenest(1-nbdy,1-nbdy,lslot),1)
c ---    ficem    - mnproc=1 returns field
         call readnest('nest/',dtime,'ficem   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,ficenest(1-nbdy,1-nbdy,lslot),1)
c ---    hsnwm(new)- mnproc=1 returns field
         call readnest('nest/',dtime,'hsnwm   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,hsnwnest(1-nbdy,1-nbdy,lslot),1)
c ---    uice     - mnproc=1 returns field
         call readnest('nest/',dtime,'uice    ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,uicenest(1-nbdy,1-nbdy,lslot),1)
c ---    vice     - mnproc=1 returns field
         call readnest('nest/',dtime,'vice    ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,vicenest(1-nbdy,1-nbdy,lslot),1)
#endif
#if defined (NEST_BIO)
#if defined (NOR05)
c ---    tracer - mnproc=1 returns field
         do ktr=1,ntracr
           write(trcnum,'(i2.2)') ktr
           do k=1,kdm
             call readnest('nest/',dtime,'trcr'//trcnum//'  ',k,
     &           fld,itdm,jtdm,imem)
             call xcaput(fld,trnest(1-nbdy,1-nbdy,k,lslot,ktr),1)
           end do
         end do
#endif /* NOR05 */
#endif /* NEST_BIO */
      end if ! baroclinic/3d fields
c
c --- Read barotropic fields (and ice)
      if (btrpflag) then 
c ---    pbavg    - mnproc=1 returns field
         call readnest('nest/',dtime,'pbavg   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,pbnest(1-nbdy,1-nbdy,lslot),1)
c ---    ubavg    - mnproc=1 returns field
         call readnest('nest/',dtime,'ubavg   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,ubnest(1-nbdy,1-nbdy,lslot),1)
c ---    vbavg    - mnproc=1 returns field
         call readnest('nest/',dtime,'vbavg   ',0, fld,itdm,jtdm,imem)
         call xcaput(fld,vbnest(1-nbdy,1-nbdy,lslot),1)
c
c ---    barotropic u-velocity at nestwalls on p-grid
!$OMP PARALLEL DO PRIVATE(j,i,hqpbot,hqpbotu,hqpbotv)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if     (ip(i,j).eq.1) then
C           hqpbot = 0.5/pbot(i,j)
C           ubpnst(i,j,lslot) = (ubnest(i,  j,lslot)*depthu(i,  j)+
C    &                           ubnest(i+1,j,lslot)*depthu(i+1,j) )
C    &                          *hqpbot
C           vbpnst(i,j,lslot) = (vbnest(i,j,  lslot)*depthv(i,j  )+
C    &                           vbnest(i,j+1,lslot)*depthv(i,j+1) )
C    &                          *hqpbot
            hqpbotu =  1./
     &         (max(1,iu(i,j)+iu(i+1,j  )) * pbot(i,j) )
            hqpbotv = 1./
     &         (max(1,iv(i,j)+iv(i  ,j+1)) * pbot(i,j) )
            ubpnst(i,j,lslot) = 
     &         (iu(i  ,j)*ubnest(i  ,j,lslot)*depthu(i  ,j)+
     &          iu(i+1,j)*ubnest(i+1,j,lslot)*depthu(i+1,j) )
     &         *hqpbotu
            vbpnst(i,j,lslot) = 
     &         (iv(i,j  )*vbnest(i,j,  lslot)*depthv(i,j  )+
     &          iv(i,j+1)*vbnest(i,j+1,lslot)*depthv(i,j+1) )
     &         *hqpbotv
          else
            pbnest(i,j,lslot) = 0.0
            ubpnst(i,j,lslot) = 0.0
            vbpnst(i,j,lslot) = 0.0
          endif
        enddo
      enddo
      end if ! btrpflag
c
c --- Rotate nest fields onto model grid directions. Outer model will 
c --- already have rotated them to east/north directions
c 
c --- baroclinic velocities
c --- TODO: xctilr before rotation
      if (bclnflag) then
         call xctilr(unest(1-nbdy,1-nbdy,1,lslot),1,kk, 1,1, halo_uv)
         call xctilr(vnest(1-nbdy,1-nbdy,1,lslot),1,kk, 1,1, halo_vv)
         do k=1,kdm
            call rotate2(unest(1-nbdy,1-nbdy,k,lslot),
     &                  vnest(1-nbdy,1-nbdy,k,lslot),
     &                  plat,plon,'l2m',.true.)
         end do
#if defined (ICE_NEST) && defined (ICE)
c ---    ice velocities
         call xctilr(uicenest(1-nbdy,1-nbdy,lslot),1,1,1,1, halo_uv)
         call xctilr(vicenest(1-nbdy,1-nbdy,lslot),1,1,1,1, halo_vv)
         call rotate2(uicenest(1-nbdy,1-nbdy,lslot),
     &               vicenest(1-nbdy,1-nbdy,lslot),
     &               plat,plon,'l2m',.true.)
#endif
      end if ! bclnflag
c --- barotropic velocities
      if (btrpflag) then
         call xctilr(ubnest(1-nbdy,1-nbdy,lslot),1,1, 1,1, halo_uv)
         call xctilr(vbnest(1-nbdy,1-nbdy,lslot),1,1, 1,1, halo_vv)
         call rotate2(ubnest(1-nbdy,1-nbdy,lslot),
     &               vbnest(1-nbdy,1-nbdy,lslot),
     &               plat,plon,'l2m',.true.)
      end if ! btrpflag
c
c --- Diagnostics option - off by default
      if (diag.and.bclnflag.and.imem==1) then
         write(ctt,'(f10.2)') dtime
         ncfil='nestfld'//trim(adjustl(ctt))//'.nc'
         call xcaget(fld,depths,0)
         call ncdraft(trim(ncfil),fld,'depths','clobber')
         call xcaget(fld,plon,0)
         call ncdraft(trim(ncfil),fld,'lon','')
         call xcaget(fld,plat,0)
         call ncdraft(trim(ncfil),fld,'lat','')
         call xcaget(fld,rmunp,0)
         call ncdraft(trim(ncfil),fld,'rmunp','')
         call xcaget(fld,rmunv,0)
         call ncdraft(trim(ncfil),fld,'rmunv','')
#if defined (ICE_NEST) && defined (ICE)
         call xcaget(fld,hicenest,0)
         call ncdraft(trim(ncfil),fld,'hicem','')
         call xcaget(fld,ficenest,0)
         call ncdraft(trim(ncfil),fld,'ficem','')
         call xcaget(fld,uicenest,0)
         call ncdraft(trim(ncfil),fld,'uicenest','')
         call xcaget(fld,vicenest,0)
         call ncdraft(trim(ncfil),fld,'vicenest','')
         call xcaget(fld,hsnwnest,0)
         call ncdraft(trim(ncfil),fld,'hsnwnest','')
#endif
#if defined (NEST_BIO)
#if defined (NOR05)
         do ktr=1,ntracr
           write(trcnum,'(i2.2)') ktr
           do k=1,kdm
             write(ckk,'(i3.3)') k
             call xcaget(fld,trnest(1-nbdy,1-nbdy,k,lslot,ktr),0)
             call ncdraft(trim(ncfil),fld,'trnest'//trcnum//'_'//ckk,'')
           end do
         end do
#endif /* NOR05 */
#endif /* NEST_BIO */
         if (bclnflag) then
         do k=1,kdm
            write(ckk,'(i3.3)') k
            call xcaget(fld,pnest(1-nbdy,1-nbdy,k,lslot),0)
            call ncdraft(trim(ncfil),fld,'pnest'//ckk,'')
            call xcaget(fld,tnest(1-nbdy,1-nbdy,k,lslot),0)
            call ncdraft(trim(ncfil),fld,'tnest'//ckk,'')
            call xcaget(fld,unest(1-nbdy,1-nbdy,k,lslot),0)
            call ncdraft(trim(ncfil),fld,'unest'//ckk,'')
         end do
         end if
         if (btrpflag) then
         call xcaget(fld,ubnest(1-nbdy,1-nbdy,lslot),0)
         call ncdraft(trim(ncfil),fld,'ubnest','')
         call xcaget(fld,vbnest(1-nbdy,1-nbdy,lslot),0)
         call ncdraft(trim(ncfil),fld,'vbnest','')
         call xcaget(fld,pbnest(1-nbdy,1-nbdy,lslot),0)
         call ncdraft(trim(ncfil),fld,'pbnest','')
         end if
      end if
      end subroutine nestcondread

c
c --- Name of new nesting  files
      function filenesting2(dtime,nestname)
      implicit none
      character(len=200) filenesting2
      real*8, intent(in) :: dtime
      character(len=*), intent(in) :: nestname
      character(len=4) :: cyy
      character(len=3) :: cdd
      character(len=2) :: chh
      integer l,trimlen, iyear, iday, ihour
      include 'common_blocks.h'
c
      call forday(dtime, yrflag, iyear,iday,ihour)
      write(cyy,'(i4.4)') iyear
      write(cdd,'(i3.3)') iday-1
      write(chh,'(i2.2)') ihour
c
      filenesting2='nest_'//cyy//'_'//cdd//'_'//chh
      if (len_trim(nestname)/=0) then
         filenesting2=trim(nestname)//trim(filenesting2)
      end if
      end function filenesting2


c --- Write one "line" of nesting data into the nest output files
c --- Files written to are header, i1, ii, j1 and jj files
      subroutine writenest(nestbase,dtime,cvar,k,fld,nx,ny,irec,imember)
      implicit none
      integer         , intent(in) :: nx, ny, irec
      character(len=*), intent(in) :: nestbase
      character(len=8), intent(in) :: cvar
      integer         , intent(in) :: k
      real,             intent(in) :: fld(nx,ny)
      real*8          , intent(in) :: dtime    
      integer         , intent(in) :: imember
c
      integer :: i1,i2,j1,j2, ibnd,j
      character(len=2) :: cbnd
      character(len=200) :: cfile
      character(len=3) :: cmem
      real*4 :: inner_ior4(nx, ny), hmin,hmax
      include 'common_blocks.h'
c
      cfile=filenesting2(dtime+1e-9,trim(nestbase))
      write(cmem,'(i3.3)') imember
c
c --- Process 4 edges
      hmin=hugel; hmax=-hugel
      do ibnd=1,4
         if     (ibnd==1)  then
            cbnd='i1'
            i1=1    ; i2=inest
            j1=1    ; j2=ny
         elseif (ibnd==2)  then
            cbnd='j1'
            i1=1    ; i2=nx
            j1=1    ; j2=inest
         elseif (ibnd==3)  then
            cbnd='ii'
            i1=nx-inest+1   ; i2=nx
            j1=1    ; j2=ny
         elseif (ibnd==4)  then
            cbnd='jj'
            i1=1    ; i2=nx
            j1=ny-inest+1   ; j2=ny
         end if
c
c ---    write to data files
         if (mnproc==1) then
            inquire(iolength=j)inner_ior4(i1:i2,j1:j2)
            inner_ior4(i1:i2,j1:j2)=fld(i1:i2,j1:j2)
            hmin=min(minval(inner_ior4),hmin)
            hmax=max(maxval(inner_ior4),hmax)
            if (irec==1) then
              if (imember==1) then
                 open(11,file=trim(cfile)//'_'//cbnd,form='unformatted',
     &             access='direct',recl=j,status='replace')
              else
                 open(11,file=trim(cfile)//'_'//cmem//'_'//cbnd,
     &           form='unformatted',access='direct',
     &           recl=j,status='replace')
              end if
            else
              if (imember==1) then
                 open(11,file=trim(cfile)//'_'//cbnd,form='unformatted',
     &             access='direct',recl=j,status='old')
              else
                 open(11,file=trim(cfile)//'_'//cmem//'_'//cbnd,
     &           form='unformatted',access='direct',recl=j,status='old')
              end if
            end if
            write(11,rec=irec) inner_ior4(i1:i2,j1:j2)
            close(11)
         end if
      end do
c --- write to header file
      if (mnproc==1 .and. imember==1) then
         if (irec==1) then
            open(10,file=trim(cfile)//'.hdr',form='formatted', 
     &              access='sequential',status='replace')
            write(10,104) dtime,kdm,inest
         else
            open(10,file=trim(cfile)//'.hdr',form='formatted', 
     &               access='sequential',status='old',
     &               position='append')
         end if
         write (10,103) cvar, k, hmin, hmax
         close(10)
      end if
  103 format(a8,": level=",i4," min/max:",2e14.6)
  104 format("dtime, kdm, inest :",f16.4,i3,i3) 
      end subroutine

c --- Read one "line" of nest data into temporary array field.
c --- Field contains data for all four edges
      subroutine readnest(nestbase,dtime,cvar,k,fld,nx,ny,imember)
      use mod_za
      implicit none
      integer         , intent(in)  :: nx, ny
      character(len=*), intent(in)  :: nestbase
      character(len=8), intent(in)  :: cvar
      integer         , intent(in)  :: k
      real            , intent(out) :: fld(nx,ny)
      real*8          , intent(in)  :: dtime           
      integer         , intent(in)  :: imember
c
      integer :: i1,i2,j1,j2, ibnd, j, tmp_kdm, tmp_inest, k2, i,
     &   ios,irec
      real*4 :: inner_ior4(nx, ny), hmin,hmax
      real*8 :: tmp_dtime
      character(len=2) :: cbnd
      character(len=8) :: cvarf
      character(len=80):: cline
      character(len=200) :: cfile
      character(len=3) :: cmem
      include 'common_blocks.h'
c
c --- Open header
      cfile=filenesting2(dtime+1e-9,trim(nestbase))
      write(cmem,'(i3.3)') imember
      if (mnproc==1) 
     &   open(10,file=trim(cfile)//'.hdr',form='formatted', 
     &           access='sequential',status='old')
      if (mnproc==1) print*, 'filename =', trim(cfile),dtime+1e-9
      call zagetc(cline,ios, 10)
      read(cline,104) tmp_dtime, tmp_kdm,tmp_inest
      if (inest/=tmp_inest .or. kdm/=tmp_kdm ) then
         if (mnproc==1) then
            write(lp,'(a)') 'Dimension error when reading nest file'
         endif
         call xcstop('(readnest)')
         stop '(readnest)'
      endif
c --- Get sought-for variable
      ios=0 ; irec=0 ; i=1
      do while (ios==0 .and. irec==0)
         !read (10,103,iostat=ios) cvarf, k2, hmin,hmax
         call zagetc(cline,ios, 10)
         if (ios .eq. 0) then
           read (cline,103,iostat=ios) cvarf, k2, hmin,hmax
         end if
         if (trim(cvar)==trim(cvarf) .and. k==k2) irec=i
         i=i+1
      end do
      if (mnproc==1) close(10)
      if (ios/=0 .or. irec==0) then
         if (mnproc==1) then
            write(lp,'(a)') 'Could not find variable '//
     &         cvar,k
         endif
         call xcstop('(readnest)')
         stop '(readnest)'
      endif
c
c --- Read 4 edges
      do ibnd=1,4
         if     (ibnd==1)  then
            cbnd='i1'
            i1=1    ; i2=inest
            j1=1    ; j2=ny
         elseif (ibnd==2)  then
            cbnd='j1'
            i1=1    ; i2=nx
            j1=1    ; j2=inest
         elseif (ibnd==3)  then
            cbnd='ii'
            i1=nx-inest+1   ; i2=nx
            j1=1    ; j2=ny
         elseif (ibnd==4)  then
            cbnd='jj'
            i1=1    ; i2=nx
            j1=ny-inest+1   ; j2=ny
         end if
c
         if (mnproc==1) then
           cfile=filenesting2(dtime+1e-9,trim(nestbase))
           inquire(iolength=j)inner_ior4(i1:i2,j1:j2)
           if (imember==1) then
              open(11,file=trim(cfile)//'_'//cbnd,form='unformatted',
     &             access='direct',recl=j,status='old')
              read(11,rec=irec) inner_ior4(i1:i2,j1:j2)
              fld(i1:i2,j1:j2)=inner_ior4(i1:i2,j1:j2)
              close(11)
           else
              open(11,file=trim(cfile)//'_'//cmem//'_'//cbnd,form=
     &        'unformatted',access='direct',recl=j,status='old')
              read(11,rec=irec) inner_ior4(i1:i2,j1:j2)
              fld(i1:i2,j1:j2)=inner_ior4(i1:i2,j1:j2)
              close(11)
           end if
         end if
      end do
c --- Diag
      if (mnproc==1) then
         write(lp,'(a,a,i3)') 'readnest: ',cvar,k
      end if
c
  103 format(a8,": level=",i4," min/max:",2e14.6)
  104 format("dtime, kdm, inest :",f16.4,i3,i3) 
      end subroutine



c --- Interpolate data in tmp arrays onto positions indicated by 
c --- inner positions (inner_...). Inner positions and bathymetri
c --- must be allocated before running this routine
      subroutine interpnest(inner_tmp,inner_tmpa,inx,iny,fld,dp,layerfl)
      use mod_xc
      implicit none
      integer, intent(in) :: inx, iny
      logical, intent(in) :: layerfl
      real, intent(out),   dimension(inx,iny) :: inner_tmp
      real, intent(inout), dimension(inx,iny) :: inner_tmpa
      real, intent(in), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   fld, dp
      real, parameter :: onem=9806. ! for now
c
      real :: ba1, ba2, ba3, ba4, dp1, dp2, dp3, dp4, dpsum,
     &        basum, rdummy
      integer :: i,j, ipiv, jpiv
      logical :: useabove, lerr
c
      lerr=.false.
#define NEST_DP_WEIGHT
c --- Do actual interpolation using what we calculated earlier 
!$OMP PARALLEL DO PRIVATE(i,j,lon,lat,ipiv,jpiv, 
!$OMP&                    ba1,ba2,ba3,ba4,
#if defined (NEST_DP_WEIGHT)
!$OMP&                    ,basum, dpsum,useabove,dp1,dp2,dp3,dp4
!$OMP&                    ) SHARED(lerr
#endif /* NEST_DP_WEIGHT */
!$OMP&                    ) 
      do j=1,inner_jdm
      do i=1,inner_idm
      inner_tmp(i,j)=0.
      if (inner_tile_mask(i,j)) then 
         ipiv=inner_ipiv(i,j)-i0 ! On local grid
         jpiv=inner_jpiv(i,j)-j0 ! On local grid
         ba1 =inner_a1  (i,j)
         ba2 =inner_a2  (i,j)
         ba3 =inner_a3  (i,j)
         ba4 =inner_a4  (i,j)

#if defined (NEST_DP_WEIGHT)
c ---    KAL -- Logic for weighting the values with layer thickness
c ---    These are the layer thicknesses surrounding this cell
         if (layerfl) then
            dp1=dp(ipiv  ,jpiv  )
            dp2=dp(ipiv+1,jpiv  )
            dp3=dp(ipiv+1,jpiv+1)
            dp4=dp(ipiv  ,jpiv+1)
         else  ! This if its a 2D var or dp
            dp1=1.1*onem
            dp2=1.1*onem
            dp3=1.1*onem
            dp4=1.1*onem
         end if
c
c ---    Summed thickness of layers
         dpsum=dp1+dp2+dp3+dp4
         useabove=.false.
c
c ---    If all grid cells are sufficiently thick, we use the
c ---    existing weights 
         if (all((/dp1,dp2,dp3,dp4/)>onem)) then
            continue
c ---    Check if averaging cells are thick enough
         elseif (dpsum>onem) then
            ! Layer thickness ok, redefine weights
            ba1=ba1*dp1
            ba2=ba2*dp2
            ba3=ba3*dp3
            ba4=ba4*dp4
            basum=ba1+ba2+ba3+ba4
c
c ---       Make sure max weight exceeds this ! value
            if ( basum<0.01*dpsum ) then
               if (layerfl) then  ! use values from layer above
                  useabove=.true.
               else
                  ! Shouldnt happen
                  lerr    =.true.
                  print *,'nestcondsave:  sum(dp*weights) too low'
               end if
            end if
c
c ---       Ensure they sum up to 1
            ba1=ba1/basum
            ba2=ba2/basum
            ba3=ba3/basum
            ba4=ba4/basum
c ---    use values from layer above
         else if (layerfl) then 
            useabove=.true.
c ---    This shouldnt happen
         else
            print *,'nestcondsave:  sum(dp) too low '
            lerr    =.true.
         end if
#endif /* NEST_DP_WEIGHT */
C
c ---    We should have a representation of this point now...
         inner_tmp(i,j)=  ba1*fld(ipiv  ,jpiv  )     
     &                   +ba2*fld(ipiv+1,jpiv  )
     &                   +ba3*fld(ipiv+1,jpiv+1) 
     &                   +ba4*fld(ipiv  ,jpiv+1)
C
#if defined (NEST_DP_WEIGHT)
c ---    use fld above
         if (useabove) then
            inner_tmp(i,j)=inner_tmpa(i,j) 
         end if
c ---    fld above, used next time
         inner_tmpa(i,j)=inner_tmp(i,j) 
#endif /* NEST_DP_WEIGHT */
      end if
      end do
      end do
c
#if defined (NEST_DP_WEIGHT)
c --- Improvized flag synchronization
      if (lerr) then
         rdummy=1.0
      else
         rdummy=0.0
      end if
      call xcmaxr(rdummy)
      lerr=nint(rdummy)==1
c
      if (lerr) then
         print *,'Error in nesting - stopping '
         call xcstop('(m_nestcondsave)')
         stop '(m_nestcondsave)'
      end if
#endif /* NEST_DP_WEIGHT */
      end subroutine


c --- Allocates variables and reads the grid of the inner model. 
c --- (Positions and bathymetry)
      subroutine readinnergrid(nestname)
      implicit none 
      character(len=*), intent(in)  :: nestname
      real*4, allocatable :: inner_io1(:)
      real*4, allocatable :: inner_io2(:,:)
      logical :: ex
      integer :: rl, n2drec
c
CKAL  Yes, I killed nestloc
c --- Read grid size from regional.grid.b
      inquire(file=trim(nestname)//'regional.grid.b', exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &           trim(nestname)//'regional.grid.b'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         open(10,file=trim(nestname)//'regional.grid.b', status='old')
         read(10,*) inner_idm
         read(10,*) inner_jdm
         close(10)
      end if
c --- inest hardcoded for now
      inner_inest=20
      inner_iidm=inner_idm-inner_inest+1
      inner_jjdm=inner_jdm-inner_inest+1
      if (mnproc==1) then
         write(lp,'(a,i3)')  '("inest     = ",i)',inner_inest
         write(lp,'(a,i3)') '("inner_idm = ",i)', inner_idm
         write(lp,'(a,i3)') '("inner_jdm = ",i)', inner_jdm
      end if
c --- Read cell center positions from regional.grid.a
      n2drec=((inner_idm*inner_jdm+4095)/4096)*4096
      allocate(inner_io1(n2drec))
      allocate(inner_io2   (inner_idm,inner_jdm))
      allocate(inner_lon   (inner_idm,inner_jdm))
      allocate(inner_lat   (inner_idm,inner_jdm))
      allocate(inner_depths(inner_idm,inner_jdm))
      inquire(file=trim(nestname)//'regional.grid.a', exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &           trim(nestname)//'regional.grid.a'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         inquire(iolength=rl) inner_io1
         open(10,file=trim(nestname)//'regional.grid.a',
     &      form='unformatted', status='old',recl=rl,access='direct')
         read(10,rec=1) inner_io2 ; inner_lon=inner_io2
         read(10,rec=2) inner_io2 ; inner_lat=inner_io2
         close(10)
      end if
c
c ---    Read the local grid
CKAL  Yes, I killed ndepths
      inquire(file=trim(nestname)//'regional.depth.a', exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &           trim(nestname)//'regional.depth.a'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         inquire(iolength=rl) inner_io2
         open(10,file=trim(nestname)//'regional.depth.a',
     &      form='unformatted', status='old',recl=rl,access='direct')
         read(10,rec=1) inner_io2
         close(10)
         inner_depths=inner_io2
      end if
c --- Mask away the hugel values
      where (inner_depths>0.5*2.0**99) inner_depths=0.
c
c --- Extend local grid to boundary
      where (inner_depths(2,:)>0.) 
     &   inner_depths(1,:)=inner_depths(2,:)
      where (inner_depths(inner_idm-1,:)>0.) 
     &   inner_depths(inner_idm,:)=inner_depths(inner_idm-1,:)
      where (inner_depths(:,2)>0.) 
     &   inner_depths(:,1)=inner_depths(:,2)
      where (inner_depths(:,inner_jdm-1)>0.)
     &   inner_depths(:,inner_jdm)=inner_depths(:,inner_jdm-1)
      deallocate(inner_io1)
      deallocate(inner_io2)
      end subroutine 


c --- NB: ipiv,jpiv refers to points on the Global (outer) grid
c --- Assign points on inner grid to points on the global grid of this model (if possible)
c --- The pivot points are now pre-processed
      subroutine readpivots(nestname)
      use mod_xc
      use mod_hycom_nersc
      implicit none
      character(len=*), intent(in)    :: nestname
c
      real :: lon,lat, lon_n, lat_n, ba1, ba2, ba3, ba4
      integer :: i,j, ipiv, jpiv,iol4
      logical :: ex, ass
      character(len=200) :: ncfil
c     character(len=11) :: tag7g
      real*4   , dimension(:,:), allocatable :: iofld4
      integer*4, dimension(:,:), allocatable :: iofldint
c
      allocate(inner_ipiv(inner_idm,inner_jdm))
      allocate(inner_jpiv(inner_idm,inner_jdm))
      allocate(inner_a1  (inner_idm,inner_jdm))
      allocate(inner_a2  (inner_idm,inner_jdm))
      allocate(inner_a3  (inner_idm,inner_jdm))
      allocate(inner_a4  (inner_idm,inner_jdm))
      allocate(inner_tile_mask(inner_idm,inner_jdm))
      allocate(iofld4    (inner_idm,inner_jdm))
      allocate(iofldint  (inner_idm,inner_jdm))
c     write(tag7g,'(i5.5,a,i5.5)')itdm,'x',jtdm
c     inquire(file=trim(nestname)//'pivots'//trim(tag7g)//'.uf',
c    &        exist=ex)
      inquire(file=trim(nestname)//'pivots.uf',
     &        exist=ex)
      inquire(iolength=iol4) iofld4
c
c --- Read pre-processed pivot points
c --- TODO: All tasks read the same file here
      if (ex) then 
c        if (mnproc==1) 
c    &      write(lp,'(a)') 'Reading pivot points from file '//
c    &          trim(nestname)//'pivots'//trim(tag7g)//'.uf'
c        open(10,file=trim(nestname)//'pivots'//trim(tag7g)//'.uf',
c    &      form='unformatted',status='old',access='direct',
c    &      recl=iol4)
         if (mnproc==1) 
     &      write(lp,'(a)') 'Reading pivot points from file '//
     &          trim(nestname)//'pivots.uf'
         open(10,file=trim(nestname)//'pivots.uf',
     &      form='unformatted',status='old',access='direct',
     &      recl=iol4)
         read(10,rec=1) iofldint; inner_ipiv = iofldint
         read(10,rec=2) iofldint; inner_jpiv = iofldint
         read(10,rec=3) iofld4  ; inner_a1   = iofld4
         read(10,rec=4) iofld4  ; inner_a2   = iofld4
         read(10,rec=5) iofld4  ; inner_a3   = iofld4
         read(10,rec=6) iofld4  ; inner_a4   = iofld4
         close(10)
c --- File missing
      else 
         if (mnproc==1)  then
            write(lp,'(a)') 'readpivots: Could not find:'
            print *,trim(nestname)//'pivots.uf'
c           print *,trim(nestname)//'pivots'//trim(tag7g)//'.uf'
         end if
         call xcstop('(nestcondsave)')
      end if
      call xcsync(flush_lp)
C
C --- netcdf diagnostics  
c     ncfil=trim(nestname)//'pivots'//trim(tag7g)//'tst.nc'
      ncfil=trim(nestname)//'pivots_tst.nc'
      call ncdraft(trim(ncfil),inner_depths,'depths','clobber')
      call ncdraft(trim(ncfil),inner_lon,'lon','')
      call ncdraft(trim(ncfil),inner_lat,'lat','')
      call ncdraft(trim(ncfil),real(inner_ipiv),'ipiv','')
      call ncdraft(trim(ncfil),real(inner_jpiv),'jpiv','')
      call ncdraft(trim(ncfil),inner_a1,'a1','')
      call ncdraft(trim(ncfil),inner_a2,'a2','')
      call ncdraft(trim(ncfil),inner_a3,'a3','')
      call ncdraft(trim(ncfil),inner_a4,'a4','')
 
c --- Local mask for this tile/process 
      inner_tile_mask=.false.
      do j=1,inner_jdm
      do i=1,inner_idm
      if (inner_depths(i,j)>0.) then 
         ipiv=inner_ipiv(i,j)
         jpiv=inner_jpiv(i,j)
c ---    Checks if pivot points are on this tile
         inner_tile_mask(i,j)= ipiv > i0 .and. ipiv <= i0+ii .and.
     &                         jpiv > j0 .and. jpiv <= j0+jj
c ---    Checks if in nesting zone of inner model
         inner_tile_mask(i,j)=inner_tile_mask(i,j).and.
     &      min(i,inner_idm-i,j,inner_jdm-j)<=inner_inest
      end if
      end do
      end do
      end subroutine
      end module mod_nesting
