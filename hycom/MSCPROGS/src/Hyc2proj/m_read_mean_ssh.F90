module m_read_mean_ssh
contains
subroutine read_mean_ssh(mean_ssh,idm,jdm,err)
implicit none
integer, intent(in)  :: idm,jdm
real,    intent(out) :: mean_ssh(idm,jdm)
logical, intent(out) :: err

integer ios
logical lmeanssh

!Read in model SSH mean
err=.false.
inquire(file='meanssh.uf',exist=lmeanssh)
if (lmeanssh) then
   open (10,file='meanssh.uf', status='unknown',form='unformatted')
   read (10,iostat=ios) mean_ssh
   close (10)
   if (ios/=0) err=.true.
else 
   print*,'m_read_meanssh: File not found'
   err=.true.
endif
end subroutine read_mean_ssh
end module m_read_mean_ssh
