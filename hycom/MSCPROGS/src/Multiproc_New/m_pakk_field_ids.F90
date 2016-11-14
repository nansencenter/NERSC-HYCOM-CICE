module m_pakk_field_ids


contains

   subroutine pakk_field_ids(fieldin,klevel,pakk_id,char2d,char3d,n2d,n3d,maxfields)
   implicit none

   character(len=*), intent(in)     :: fieldin
   character(len=5), intent(out)    :: pakk_id
   integer,          intent(in)     :: maxfields
   integer,          intent(inout)  :: n2d,n3d,klevel
   character(len=5), intent(out), dimension(maxfields) :: char2d, char3d

   integer :: i2d,i3d
   logical :: match


   ! Registered ave field ids translated to pakk ids
   ! 2D Vars
   pakk_id='EMPTY'
   if (trim(fieldin)=='ssh') then
      pakk_id='SSH'
   elseif (trim(fieldin)=='salflx') then
      pakk_id='SALFL'
   elseif (trim(fieldin)=='surflx') then
      pakk_id='RSURF'
   elseif (trim(fieldin)=='rsalflx') then
      pakk_id='RSALF'
   elseif (trim(fieldin)=='rsurflx') then
      pakk_id='SURFL'
   elseif (trim(fieldin)=='uice') then
      pakk_id='UICE'
   elseif (trim(fieldin)=='vice') then
      pakk_id='VICE'
   elseif (trim(fieldin)=='spdice') then
      pakk_id='ISPD'

   ! Derived 2D vars
   elseif (trim(fieldin)=='UBAVG') then
      pakk_id='UBAVG'
   elseif (trim(fieldin)=='VBAVG') then
      pakk_id='VBAVG'

   ! Icestate heat fluxes over leads
   elseif (trim(fieldin)=='lead_tot') then
      pakk_id='LFTOT'
   elseif (trim(fieldin)=='lead_sw') then
      pakk_id='LFSW' 
   elseif (trim(fieldin)=='lead_lw') then
      pakk_id='LFLW' 
   elseif (trim(fieldin)=='lead_trb') then
      pakk_id='LFTRB' 

   ! Icestate heat fluxes over ice
   elseif (trim(fieldin)=='ice_grw') then
      pakk_id='IGROW'
   elseif (trim(fieldin)=='ice_ntop') then
      pakk_id='INTOP' 
   elseif (trim(fieldin)=='ice_nbot') then
      pakk_id='INBOT' 
   elseif (trim(fieldin)=='ice_swfl') then
      pakk_id='ISWFL' 



   ! Icestate derived 2D vars
   elseif (trim(fieldin)=='ficem') then
      pakk_id='FICEM'
   elseif (trim(fieldin)=='hicem') then
      pakk_id='HICEM'


   ! 3D Vars
   elseif (trim(fieldin)=='u-vel.' .or. trim(fieldin)=='utot' ) then
      pakk_id='UT'
   elseif (trim(fieldin)=='v-vel.' .or. trim(fieldin)=='vtot' ) then
      pakk_id='VT'
   elseif (trim(fieldin)=='thknss' .or. trim(fieldin)=='pres' ) then
      pakk_id='DP   '
   elseif (trim(fieldin)=='temp') then
      pakk_id='TEM'
   elseif (trim(fieldin)=='salin'  .or. trim(fieldin)=='saln') then
      pakk_id='SAL'
   elseif (trim(fieldin)=='kinetic') then
      pakk_id='SPEED'

   ! Derived 3D vars
   elseif (trim(fieldin)=='EKE') then
      pakk_id='EKE'
   elseif (trim(fieldin)=='MKE') then
      pakk_id='MKE'

   ! Icestate 3D vars
   elseif (trim(fieldin)=='hice') then
      pakk_id='HICE'
   elseif (trim(fieldin)=='fice') then
      pakk_id='FICE'
   elseif (trim(fieldin)=='hsnw') then
      pakk_id='HSNW'
   end if

   ! Return if no match in registered ave field ids
   if (pakk_id == 'EMPTY') return


   ! Check if field is in char2d  or char3d array
   match=.false.
   if (klevel==0) then

      do i2d=1,n2d
         match=match.or.trim(char2d(i2d))==trim(pakk_id)
      end do

      ! Add to char array if no match
      if (.not.match) then
         n2d=n2d+1
         char2d(n2d)=pakk_id
      end if

   elseif (klevel>0) then
      do i3d=1,n3d
         match=match.or.trim(char3d(i3d))==trim(pakk_id)
      end do

      ! Add to char array if no match
      if (.not.match) then
         n3d=n3d+1
         char3d(n3d)=pakk_id
      end if
   end if

   end subroutine
end module
      
         
      
         

   
