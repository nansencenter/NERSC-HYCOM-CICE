module m_densitylevels
contains
subroutine densitylevels(fname,maxA,nlev,cid)
   implicit none
   character(len=100),  intent(in) :: fname
   integer,             intent(in) :: nlev
   real,                intent(in) :: maxA             
   character(len=1),    intent(in) :: cid

   integer i,j
   real da
   character(len=100) outname

   print *,trim(fname)
   print *,cid
   print *,nlev
   print *,maxA

   outname=' '
   i=15
   outname(1:i)=fname(1:i)
   select case (cid)
   case ('U')
      outname(i+1:i+1+8)='_lev.mcr'
   case ('F')
      outname(i+1:i+1+10)='_F_lev.mcr'
   case ('T')
      outname(i+1:i+1+10)='_T_lev.mcr'
   end select

   print *,'+++',trim(outname),'+++'


   open(10,file=trim(outname))
      da=maxA/float(nlev-1)
      write(10,'(a)')'#!MC 700'
      write(10,'(a)')'$!MACROFUNCTION'
      select case (cid)
      case('U')
         write(10,'(a,a,a)')'  NAME = "',trim(fname),'_lev"'
      case('F')
         write(10,'(a,a,a)')'  NAME = "',trim(fname),'_F_lev"'
      case('T')
         write(10,'(a,a,a)')'  NAME = "',trim(fname),'_T_lev"'
      end select
      write(10,'(a)')'$!CONTOURLEVELS NEW'
      write(10,'(a)')'  RAWDATA'
      write(10,'(i2)')nlev
      do i=1,nlev
         write(10,'(f12.8)')0.0+da*float(i-1)
      enddo
      write(10,'(a)')'$!ENDMACROFUNCTION'
   close(10)

end subroutine densitylevels
end module m_densitylevels
