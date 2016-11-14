module m_interpug
contains
! --- -----------------------------------------------------------------
! --- Replacement of horizontal interpolation routines bilin and bicubic
! --- with a single routine.
! --- 
! --- "interpug" - interpolation routine for "uniformly spaced" (in lon
! --- lat space) and global data sets which are periodic in longitude. 
! --- Interpolates data from data grid (old) to input grid (new) with 
! --- positions specified by newlon, newlat
! --- NB: Also handles Gaussian grids from NCEP. So, not "uniform" in latitude
! ---
! --- Interpolation methods:
! --- itype == 0 : bilinear interpolation. Continous, but discontinuous
! --- 1st derivatives. Always monotonic
! ---
! --- itype == 1 : bicubic interpolation. Continous 0-th and 1st order 
! --- derivatives. 2nd order cross derivative continous at corners. Not 
! --- Monotonic
! ---
! --- Assumptions: 
! ---  1) Data grid is uniform in lon and lat directions. Grid indices
! ---     increase as longitude/latitude increases.
! ---  OR Data grid is uniform in lon direction. Grid indices increase as 
! ---     longitude increases. Latitude uses "Gaussian" points (right now
! ---     any monotonically increasing/decreasing latitude vector will
! ---     work). 
! ---  2) Grid is periodic in longitude direction. Should cover the 
! ---     globe. 
! --- -----------------------------------------------------------------
      subroutine interpug(old,onx,ony,olonref,olatref,odlon,odlat,new,newlon,newlat,itype,gausslat)
      use mod_xc
      implicit none
! --- Dims of old (data) grid, as well as reference lon/lat and grid increment
      integer, intent(in) :: onx, ony,itype
      real, intent(in)    :: olonref, olatref    
      real, intent(in)    :: odlon  , odlat
      real, intent(in)    :: old(onx,ony)! old grid
!
! --- Longitude/Latitude of new grid, as well as new (interpolated) values
      real, intent(in), dimension (1:idm,1:jdm) :: newlon, newlat
      real, intent(out)   :: new  (1:idm,1:jdm)
!
! --- Optional argument indicating gaussian grid
      real, intent(in), optional  :: gausslat(ony)
!
      integer i,j,ipos,jpos,ia,ib,ja,jb,ib2,jb2,numerr,j2
      real    aa,bb,lon,minlat,rnumerr, newlon2,dlatg
      real :: f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,fxy00,fxy10,fxy01,fxy11, a1,a2,a3,a4
      real :: rhs(16), coeffs(16), radian
      logical :: gaussgrid
      real, parameter :: thlat=30. ! TODO: tune this parameter
      radian=asin(1.)/90.
!
! --- Security check
      if (itype<0 .and. itype>1) then
         write(lp,'(a)') 'Invalid option for itype'
         call xcstop('(mod_hycom_nersc:interpug)')
      end if
!
! --- Check for gaussian
      gaussgrid=.false.
      if (present(gausslat)) then
         gaussgrid=.true.
      end if
!
! --- Start interpolation
      numerr=0
      do j=1,jdm
      do i=1,idm
! ---    New index in old data. New point is between
! ---    ipos, ib and jpos, jb
         newlon2 = mod(newlon(i,j)-olonref+360.d0,360.d0) + olonref
         ipos =int((newlon2-olonref)/odlon+1.0d0)  ! [1:onx]
         if (ipos>onx .or. ipos < 1 ) then
            write(lp,'(a,2i5,2f10.2)') 'ipos error ', &
            ipos,onx,newlon2,olonref+odlon*(onx-1)
            numerr=numerr+1
         endif
         ia =mod(onx+ipos-2,onx)+1    ! Periodic
         ib =mod(ipos,onx)+1          ! Periodic
         ib2=mod(ib  ,onx)+1          ! Periodic
!
! ---    jpos on uniform grid
         if (.not.gaussgrid) then
            jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
            if (jpos>ony .or. jpos < 1 ) then
               write(lp,'(a,2i5,2f10.2)') 'jpos error ',&
                 jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
               numerr=numerr+1
            endif
! ---       TODO: Latitude "Wrap-over" at North Pole ?
            ja   =min(max(1,jpos-1),ony)
            jb   =min(max(1,jpos+1),ony)
            jb2  =min(max(1,jpos+2),ony)
!
! ---    jpos on gaussian grid. NB: There exists formulas for the
! ---    Gaussian locations, and they should be used 
         else
            if (odlat>0.) then ! latitude inrease with increasing index
               jpos=ony
               do j2=1,ony
                  if (gausslat(j2)>newlat(i,j)) then
                     jpos=j2
                     exit
                  end if
               end do
! ---          TODO: Latitude "Wrap-over" at North Pole ?
               jpos =min(max(1,jpos  ),ony)
               ja   =min(max(1,jpos-1),ony)
               jb   =min(max(1,jpos+1),ony)
               jb2  =min(max(1,jpos+2),ony)
            else ! Latitude decreases with increasing index
               jpos=1
               do j2=1,ony
                  if (gausslat(j2)<newlat(i,j)) then
                     jpos=j2
                     exit
                  end if
               end do
! ---          TODO: Latitude "Wrap-over" at North Pole ?
               jpos =min(max(1,jpos  ),ony)
               ja   =min(max(1,jpos+1),ony)
               jb   =min(max(1,jpos-1),ony)
               jb2  =min(max(1,jpos-2),ony)
            end if
         endif
!        if (j==250) print '(a,5i4)','i',i,ia,ipos,ib,ib2
!
! ---    Grid distance new point -> ipos, jpos. aa,bb in [0,1]
         aa=(newlon2 - olonref-real(ipos-1)*odlon)/odlon
         if (.not. gaussgrid) then
            bb=(newlat(i,j) - olatref-real(jpos-1)*odlat)/odlat
         else
            if (gausslat(jb)==gausslat(jpos)) then
               dlatg=abs(90-gausslat(jpos))
            else
               dlatg=gausslat(jb)-gausslat(jpos)
            end if
            dlatg =abs(dlatg)
            bb=(newlat(i,j) - gausslat(jpos))/dlatg
         end if
!
! ---    Catch errors - but dont stop until after loop
         if ((aa > 1.0).or.(aa < 0.0)) then
            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
            print *,'interpug: invalid aa',aa
            numerr=numerr+1
         endif
         if ((bb > 1.0).or.(bb < 0.0)) then
            if (gaussgrid) then
               write(*,'(3i5,3f10.2)')i,j,jpos,gausslat(jpos),&
                  newlat(i,j),gausslat(jb)
            else
               write(*,'(3i5,3f10.2)')i,j,jpos,lon,newlat(i,j),&
                 olatref+(jpos-1)*odlat
            end if
            print *,'interpug: invalid bb',bb
            numerr=numerr+1
         endif
!
! ---   Set up bilinear weights
        if (itype==0) then
! ---      Bilinear weights
           a1=(1.0-aa)*(1.0-bb)
           a2=aa*(1.0-bb)
           a3=aa*bb
           a4=(1.0-aa)*bb
! ---      New data value
           new(i,j) = a1*old(ipos,jpos)+a2*old(ib  ,jpos)+&
              a3*old(ib  ,jb  )+a4*old(ipos,jb  )
! ---   TODO: most of this can be re-used if (ipos,jpos) hasnt changed
! ---   Set up function and derivatives at ecmwf nodes
! ---   TODO: Change plat in derivatives to data grid lat values
        elseif (itype==1) then
          f00  =old(ipos,jpos)
          f10  =old(ib  ,jpos)
          f01  =old(ipos,jb  )
          f11  =old(ib  ,jb  )
! ---     X derivative with gridspacing 1 - LB: no gridsize needed
          fx00 = 0.5*(old(ib  ,jpos) - old(ia  ,jpos))
          fx10 = 0.5*(old(ib2 ,jpos) - old(ipos,jpos))
          fx01 = 0.5*(old(ib  ,jb  ) - old(ia  ,jb  ))
          fx11 = 0.5*(old(ib2 ,jb  ) - old(ipos,jb  ))
! ---     Y derivative with gridspacing 1 
          fy00 = 0.5*(old(ipos,jb  ) - old(ipos,ja  ))
          fy10 = 0.5*(old(ib  ,jb  ) - old(ib  ,ja  ))
          fy01 = 0.5*(old(ipos,jb2 ) - old(ipos,jpos))
          fy11 = 0.5*(old(ib  ,jb2 ) - old(ib  ,jpos))
! ---     Cross derivative with gridspacing 1  
          fxy00=0.25*( old(ib  ,jb )-old(ib  ,ja  )- &
            (old(ia  ,jb )-old(ia  ,ja  )))
          fxy10=0.25*( old(ib2 ,jb )-old(ib2 ,ja  )- &
            (old(ipos,jb )-old(ipos,ja  )))
          fxy01=0.25*( old(ib  ,jb2)-old(ib  ,jpos)- &
            (old(ia  ,jb2)-old(ia  ,jpos)))
          fxy11=0.25*( old(ib2 ,jb2)-old(ib2 ,jpos)- &
            (old(ipos,jb2)-old(ipos,jpos)))
! ---     RHS of coeff equation
          rhs=(/f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11, &
             fxy00,fxy10,fxy01,fxy11/)
! ---     Solve matrix for cubic coeffs
! ---     TODO: optimize this routine
          coeffs=cubiccoeff(rhs)
! ---     Calculate solution
          new(i,j)=cubicsol(coeffs,aa,bb)
        end if
      end do
      end do
! --- Halt on errors
      if (numerr>0) then
         write(lp,'(a)')'Error(s) occured in interpug..'
         call xcstop('(interpug)')
         stop '(interpug)'
      end if
      end subroutine interpug

      function cubiccoeff(rhs)
      implicit  none
      real, dimension(16), intent(in) :: rhs
      real, dimension(16) :: cubiccoeff
! --- TODO: Matrix is sparse - so there is room for reducing the number
! --- of operations (i.e avoid matmul)
      real, dimension(16*16), parameter :: &
     invcb=(/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,&
              0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,&
              0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,&
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,&
              0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,&
              0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,&
              0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,&
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,&
              0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,&
              0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,&
              0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,&
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,&
              0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,&
              0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1,&
              0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,&
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1/)
      real,parameter,dimension(16,16):: invcb2=reshape(invcb,(/16,16/))
      cubiccoeff=matmul(invcb2,rhs)
      end function
!  --- Final calculation of bicubic interpolation

      real function cubicsol(coeffs,aa,bb)
      implicit none
      real, intent(in) :: aa,bb,coeffs(16)
      real :: coeffs2(4,4)
      integer :: i,j
! --- Reshape to c_ij format (see any work on bicubic interpolation)
      coeffs2=reshape(coeffs,(/4,4/))
      cubicsol=0.
      do j=1,4
      do i=1,4
         cubicsol=cubicsol+coeffs2(i,j)*aa**(i-1)*bb**(j-1)
      end do
      end do
      end function cubicsol
end module m_interpug
