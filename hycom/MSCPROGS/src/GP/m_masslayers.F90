module m_masslayers
contains
subroutine masslayers(z,nz,ipp,nip)
! setting up vectors to interpolate from neglecting thin layers
! For spline version
      implicit none
      integer               , intent(in)  :: nz
      real,    dimension(nz), intent(in)  :: z
      integer, dimension(nz), intent(out) :: ipp
      integer               , intent(out) :: nip
      integer l,k

      ipp(1)=1
      l=1
      do k=2,nz
         if ((z(k)-z(k-1)) > 0.05) then
            l=l+1
            ipp(l)=k
         endif
      enddo
      nip=l

end subroutine masslayers
end module m_masslayers
