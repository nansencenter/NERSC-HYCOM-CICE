module mod_curviint_interp

   integer, dimension(:,:), allocatable, save :: &
      ipiv,jpiv
   logical, dimension(:,:), allocatable, save :: &
      ass
   real   , dimension(:,:,:), allocatable, save :: &
      a
   logical, save, private :: periodic

contains


   ! Sets up interpolation details (pivot points, bilinear weights)
   subroutine interp_setup(lmsk,llon,llat,gmsk,glon,glat,lperiodic)
   use mod_xc_local , only : lnx=>idm,lny=>jdm
   use mod_xc_global, only : gnx=>idm,gny=>jdm
   implicit none
   logical, dimension(lnx,lny), intent(in) :: lmsk
   logical, dimension(gnx,gny), intent(in) :: gmsk
   real, dimension(lnx,lny), intent(in) :: llon,llat
   real, dimension(gnx,gny), intent(in) :: glon,glat

   logical, intent(in) :: lperiodic
   logical, dimension(lnx,lny) :: donep
   real   , dimension(lnx,lny) :: pass
   integer :: i,j

   periodic=lperiodic

   ! Bilinear interpolation setup
   donep=.false.
   call bilin_setup(llon,llat,lmsk,donep,     & 
                    glon,glat,gmsk)
   do j=1,lny
   do i=1,lnx
      if (donep(i,j)) then
         pass(i,j)=1
      else
         pass(i,j)=0
      end if
   end do
   end do
   print *,'bilin_confgrd done'

   ! Now remaining wet points except boundary points are filled in by extrapolation
   call extrapolate2(lmsk,donep)
   do j=1,lny
   do i=1,lnx
      if (donep(i,j).and.pass(i,j)==0) then
         pass(i,j)=2
      end if
   end do
   end do
   write(*,*)'extrapolate done'



   ! Check if any points are not properly set up
   if (any((.not.donep).and.lmsk)) then
      print '(a,i6,a)','failed to set up ',count((.not.donep).and.lmsk), &
      ' points .. Check for isolated points on the bathymetry grid '
      print '(a)','See error messages above for candidate points'
      stop '(curviint)'
   end if
   end subroutine



   ! Bilinear calc - uses bilinear interpolation, and it also
   ! uses layer thickness when interpolatiing
   subroutine bilinear_calc(cfld,lfld,lofld,lmsk,gfld,dp,lcoord)
   use mod_xc_local , only : lnx=>idm,lny=>jdm
   use mod_xc_global, only : gnx=>idm,gny=>jdm
   implicit none
   character(len=*), intent(in) :: cfld
   real, dimension(lnx,lny), intent(in) :: lofld
   real, dimension(lnx,lny), intent(out) :: lfld
   logical, dimension(lnx,lny), intent(in) :: lmsk
   real, dimension(gnx,gny), intent(in) :: gfld
   real, dimension(gnx,gny), intent(in) :: dp
   integer, intent(in) :: lcoord

   integer :: i,j,i2,j2,i2pb,j2pb
   real    :: ba1, ba2, ba3, ba4, basum
   real, parameter :: onem=9806.
   integer :: dbgout=-1


   ! Interpolate and  extrapolate ocean variables
   do j=1,lny
   do i=1,lnx
   if (lmsk(i,j)) then
      i2=ipiv(i,j)
      j2=jpiv(i,j)
      if (periodic) then
         i2pb=mod(i2,gnx)+1
      else
         i2pb=i2+1
      end if
      j2pb=j2+1

      if (i2pb/=i2+1) then
         print *,'NB:periodic grid set up but needs testing'
         stop
      end if
      
      ! Straight-forward bilinear interpolation for these fields
      if (trim(cfld)/='temp' .and. trim(cfld)/='saln' .and.  trim(cfld)/='th3d') then
         lfld(i,j)=                     &
            gfld(i2  ,j2  )*a(i,j,1) +  &
            gfld(i2pb,j2  )*a(i,j,2) +  &
            gfld(i2pb,j2pb)*a(i,j,3) +  &
            gfld(i2  ,j2pb)*a(i,j,4) 
         if (trim(cfld)=='dp' .and. lfld(i,j)<0.) then
            if (dbgout<0) then
               dbgout=99
               open(unit=dbgout,file='debug.out',status='replace',form='formatted')
            else
               write(dbgout,*) '***Got negative dp ... ',i,j,lcoord
               write(dbgout,*) lfld(i,j)
               write(dbgout,*) a(i,j,:)
               write(dbgout,*) gfld(i2,j2)
               write(dbgout,*) gfld(i2pb,j2)
               write(dbgout,*) gfld(i2pb,j2pb)
               write(dbgout,*) gfld(i2,j2pb)
            end if
            lfld(i,j)=0.
         end if
      else ! temp, salt and density needs to be layer weighted as well

         ba1=a(i,j,1)*dp(i2  ,j2  )
         ba2=a(i,j,2)*dp(i2pb,j2  )
         ba3=a(i,j,3)*dp(i2pb,j2pb)
         ba4=a(i,j,4)*dp(i2  ,j2pb)
         basum=ba1+ba2+ba3+ba4
         if (basum>=onem) then ! Layer weighted value
            lfld(i,j)=                     &
               gfld(i2  ,j2  )*ba1/basum +  &
               gfld(i2pb,j2  )*ba2/basum +  &
               gfld(i2pb,j2pb)*ba3/basum +  &
               gfld(i2  ,j2pb)*ba4/basum 
         else if (basum<onem.and.lcoord>1) then ! Use layer value above
            lfld(i,j)=lofld(i,j)
         else if (basum<onem) then ! abort
            print *,'basum error !'
            print *,i,j,lcoord
            print *,i2,i2pb,j2,j2pb
            print *,'****'
            print *,a(i,j,1),dp(i2  ,j2  )
            print *,a(i,j,2),dp(i2pb,j2  )
            print *,a(i,j,3),dp(i2pb,j2pb)
            print *,a(i,j,4),dp(i2  ,j2pb)
            !
            stop '(curviint)'
         end if
      end if
   end if
   end do
   end do
   end subroutine



! private routine for setting up bilinear weights. 
subroutine bilin_setup(llon,llat,lmsk,donep,    &
                         glon,glat,gmsk)
   use mod_xc_local , only : lnx=>idm,lny=>jdm
   use mod_xc_global, only : gnx=>idm,gny=>jdm
   use mod_confmap
   implicit none
   real,     intent(in)    :: llon(lnx,lny)
   real,     intent(in)    :: llat(lnx,lny)
   logical,  intent(in)    :: lmsk(lnx,lny)
   logical,  intent(out)   :: donep(lnx,lny)
   real,     intent(in)    :: glon(gnx,gny)
   real,     intent(in)    :: glat(gnx,gny)
   logical,  intent(in)    :: gmsk(gnx,gny)

   integer l,i,j,ia,ib,k
   integer ipib,jpib
   real aa,bb,atmp(4)
   logical masktmp(4)
   integer :: ip2,jp2
   real lat_n,lon_n
   real dpsum

   real, parameter:: epsil=1.0e-11
   real t1,t2,s1,s2,r1,r2
   real, parameter ::  onem=9806.
   real, parameter ::  onemm=98.06
   integer iloc(2), igrace, jgrace
   logical inirange

   allocate(ipiv    (lnx,lny))   ! i pivot point
   allocate(jpiv    (lnx,lny))   ! j pivot point
   allocate(ass     (lnx,lny))   ! Flag
   allocate(a       (lnx,lny,4)) ! bilinear weights


   ! For each gridpoint in local grid if wet
   do j=1,lny
   do i=1,lnx
      if (lmsk(i,j)) then 
         t1=llon(i,j)
         t2=llat(i,j)
         !print *,i,j,t1,t2
         call oldtonew(t2,t1,lat_n,lon_n)
         !call newtoold(lat_n,lon_n,s2,s1)
         call pivotp(lon_n,lat_n,ipiv(i,j),jpiv(i,j))

         if (periodic) then
            ipib=mod(ipiv(i,j),gnx)+1
            inirange=.true.
         else
            ipib=ipiv(i,j)+1
            inirange=ipiv(i,j)>=1 .and. ipiv(i,j) < gnx
         end if
         jpib=jpiv(i,j)+1

         ! grace for i
         igrace=min(gnx-ipiv(i,j),ipiv(i,j)-1) ! negative when ipiv < 1 or ipiv > nxl
         ! grace for j
         jgrace=min(gny-jpiv(i,j),jpiv(i,j)-1) ! negative when jpiv < 1 or jpiv > nxl

         !! Brutally stop if pivot points are outside of grid
         !if (ipiv(i,j)>gnx-1 .or. ipiv(i,j)<1 .or. &
         !    jpiv(i,j)>gny-1 .or. jpiv(i,j)<1 ) then
         !   print *,'pivot point is outside global grid!'
         !   print *,'pivot point  (i,j) :',ipiv(i,j),jpiv(i,j)
         !   print *,'global grid dim    :',gnx,gny
         !   stop '(bilin_confgrid)'
         !end if




         if (inirange .and. jpiv(i,j)>=1 .and. jpiv(i,j)<gny) then

            if (gmsk(ipiv(i,j),jpiv(i,j)) .and. gmsk(ipib,jpiv(i,j)) .and.&
                gmsk(ipiv(i,j),jpib     ) .and. gmsk(ipib,jpib    )) then
               call bilincoeff(glon,glat,gnx,gny,llon(i,j),llat(i,j), &
                               ipiv(i,j),jpiv(i,j),atmp(1),atmp(2),atmp(3),atmp(4))
               a(i,j,:)=atmp
               ass(i,j)=.true.
            else
               ass(i,j)=.false. 
               a(i,j,:)=0
               ipiv(i,j)=0; jpiv(i,j)=0
            endif

         else if ( igrace> -3 .and. jgrace>-3) then
            print '(a,2i5,a,2i5,a)','Point ',i,j,'with pivot ',ipiv(i,j),jpiv(i,j),' saved by grace '
            ass(i,j)=.false. 
            a(i,j,:)=0
            ipiv(i,j)=0; jpiv(i,j)=0
         else
            print '(a,2i5,a,2i5,a)','Point ',i,j,'with pivot ',ipiv(i,j),jpiv(i,j),' not saved by grace '
            stop '(bilin_confgrd)'
         endif
      end if
      if (ass(i,j)) then 
         donep(i,j)=.true.
      endif
   enddo
   enddo
end subroutine bilin_setup


! private routine for extrapolating data - sets up pivot points and weights
subroutine extrapolate2(lmsk,donep)
use mod_xc_local , only : lnx=>idm,lny=>jdm
implicit none
logical, intent(in)    :: lmsk   (lnx,lny)
logical, intent(inout) :: donep  (lnx,lny)
integer l,i,j,ia,ib,isign,ifalse,jja,iia,index,k,ja
logical donetmp(lnx,lny)
! Now remaining wet points except bounday points are filled in by extrapolation

   do l=1,200
      if (mod(l,2) == 0) isign=1
      if (mod(l,2) == 1) isign=-1

      ifalse=0
      donetmp=.false.
      do j=1,lny
      do i=1,lnx
         if (lmsk(i,j).and.(.not.donep(i,j))) then   
            ifalse=ifalse+1
            do jja=-1,1
               ja=j+jja*isign

               !if (ja < 1) ja=lny
               !if (ja > lny) ja=1
               if (ja < 1) ja=1
               if (ja > lny) ja=lny

               do iia=-1,1
                  ia=i+iia*isign

                  if (periodic) then
                     ia=mod(lnx+ia-1,lnx)+1
                  else
                     if (ia > lnx) ia=lnx
                     if (ia < 1) ia=1
                  end if

                  !if (i==11.and.j==96) print *,ia,ja

                  if (lmsk(ia,ja).and.donep(ia,ja)) then
                     ipiv(i,j)  =ipiv(ia,ja)
                     jpiv(i,j)  =jpiv(ia,ja)
                     a   (i,j,:)= a (ia,ja,:)
                     donetmp(i,j)=.true.
                     ifalse=ifalse-1
                     exit
                  endif
               enddo
               if (donetmp(i,j)) exit
            enddo
         endif
      enddo
      enddo
      print *,'ifalse after iteration',l,ifalse
      where (donetmp) donep=.true.
      if (ifalse == 0) exit
   enddo
   
   if (ifalse /= 0) then
      do j=1,lny
      do i=1,lnx
         if (lmsk(i,j).and.(.not.donep(i,j))) then
            print '(a,2I6)',' Ifalse problem in grip point (i,j): ',i,j
         endif
      enddo
      enddo
   endif
end subroutine extrapolate2



end module 
