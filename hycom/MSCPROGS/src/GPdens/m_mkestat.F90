module m_mkestat
contains

subroutine mkestat(cm,gp,cmfilt,gpfilt,cmtide,gptide,nrm,deep)
   use mod_data_new
   use mod_netcdf_ops
   implicit none
   integer, intent(in) ::  nrm
   integer, intent(in) ::  deep
   type(data_new), intent(in) :: cm(nrm),gp(nrm)
   type(data_new), intent(in) :: cmfilt(nrm),gpfilt(nrm)
   type(data_new), intent(in) :: cmtide(nrm),gptide(nrm)

   type(data_new) cmmeso(nrm),gpmeso(nrm)


   logical ex
   integer m
   real cm_mketot,cm_mkefull,cm_mketide,cm_mkemeso,cm_mkemean,cmmean_u,cmmean_v
   real gp_mketot,gp_mkefull,gp_mketide,gp_mkemeso,gp_mkemean,gpmean_u,gpmean_v

   real runmeanGP(nrm),runtideGP(nrm),runmesoGP(nrm),runfullGP(nrm)
   real runmeanCM(nrm),runtideCM(nrm),runmesoCM(nrm),runfullCM(nrm)
   integer :: ncid, vid, vertind


   cmmean_u=0.0
   cmmean_v=0.0
   gpmean_u=0.0
   gpmean_v=0.0
   do m=1,nrm
      cmmean_u=cmmean_u+cm(m)%u
      cmmean_v=cmmean_v+cm(m)%v
      gpmean_u=gpmean_u+gp(m)%u
      gpmean_v=gpmean_v+gp(m)%v
      runmeanCM(m)=sqrt( (cmmean_u/float(m))**2 + (cmmean_v/float(m))**2 )   
      runmeanGP(m)=sqrt( (gpmean_u/float(m))**2 + (gpmean_v/float(m))**2 )   
   enddo
   cmmean_u=cmmean_u/float(nrm)
   cmmean_v=cmmean_v/float(nrm)
   gpmean_u=gpmean_u/float(nrm)
   gpmean_v=gpmean_v/float(nrm)

   cm_mkemean=sqrt(cmmean_u**2+cmmean_v**2)
   gp_mkemean=sqrt(gpmean_u**2+gpmean_v**2)


   do m=1,nrm
      cmmeso(m)%u=cmfilt(m)%u-cmmean_u
      cmmeso(m)%v=cmfilt(m)%v-cmmean_v
      gpmeso(m)%u=gpfilt(m)%u-gpmean_u
      gpmeso(m)%v=gpfilt(m)%v-gpmean_v
   enddo

   cm_mkefull=0.0
   gp_mkefull=0.0
   cm_mketide=0.0
   gp_mketide=0.0
   cm_mkemeso=0.0
   gp_mkemeso=0.0
   do m=1,nrm
      cm_mkefull=cm_mkefull+cm(m)%u**2+cm(m)%v**2
      gp_mkefull=gp_mkefull+gp(m)%u**2+gp(m)%v**2

      cm_mketide=cm_mketide+cmtide(m)%u**2+cmtide(m)%v**2
      gp_mketide=gp_mketide+gptide(m)%u**2+gptide(m)%v**2

      cm_mkemeso=cm_mkemeso+cmmeso(m)%u**2+cmmeso(m)%v**2
      gp_mkemeso=gp_mkemeso+gpmeso(m)%u**2+gpmeso(m)%v**2

      runmesoCM(m)=sqrt( cm_mkemeso/float(m) )   
      runmesoGP(m)=sqrt( gp_mkemeso/float(m) )   

      runtideCM(m)=sqrt( cm_mketide/float(m) )   
      runtideGP(m)=sqrt( gp_mketide/float(m) )   

      runfullCM(m)=sqrt( cm_mkefull/float(m) )   
      runfullGP(m)=sqrt( gp_mkefull/float(m) )   
   enddo
   cm_mkefull=sqrt(cm_mkefull/float(nrm))
   gp_mkefull=sqrt(gp_mkefull/float(nrm))
   cm_mketide=sqrt(cm_mketide/float(nrm))
   gp_mketide=sqrt(gp_mketide/float(nrm))
   cm_mkemeso=sqrt(cm_mkemeso/float(nrm))
   gp_mkemeso=sqrt(gp_mkemeso/float(nrm))

   cm_mketot=sqrt(cm_mketide**2+cm_mkemeso**2+cm_mkemean**2)
   gp_mketot=sqrt(gp_mketide**2+gp_mkemeso**2+gp_mkemean**2)



!   inquire(file='mkestat.dat',exist=ex)
!   if (.not.ex) then
!      open(10,file='mkestat.dat')
!         write(10,'(a,i3,a)')'{\\bf Mooring ',station,&
!           &'} #  \\multicolumn{5}{c}{Model} # \\multicolumn{5}{c}{Measured} \\\\'
!
!         write(10,'(16(a,1x))')'Depth #','KE(full) #','KE(tot)  #','KE(tide) #','KE(meso) #','KE(mean) #',&
!                                        &'KE(full) #','KE(tot)  #','KE(tide) #','KE(meso) #','KE(mean) \\\\'
!      close(10)
!   endif
!
!
!   open(10,file='mkestat.dat',position='append')
!      write(10,'(i5,1x,10(a,f10.2,1x),a)')deep,&
!        '#',gp_mkefull,'#',gp_mketot,'#',gp_mketide,'#',gp_mkemeso,'#',gp_mkemean,&
!        '#',cm_mkefull,'#',cm_mketot,'#',cm_mketide,'#',cm_mkemeso,'#',cm_mkemean,'\\\\'
!   close(10)
!
!   open(10,file='runmean.dat')
!      do m=1,nrm
!         write(10,'(7f10.2)')float(m)/24.0,runtideGP(m),runtideCM(m),runmesoGP(m),runmesoCM(m),runmeanGP(m),runmeanCM(m)
!      enddo
!   close(10)


   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)
   call ncinqdefvertvar(ncid,vid,deep,vertind)
   call ncinqputvar(ncid,'mke_full_GP',(/vid/),(/gp_mkefull/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_tot_GP' ,(/vid/),(/gp_mketot /),start=(/vertind/))
   call ncinqputvar(ncid,'mke_tide_GP',(/vid/),(/gp_mketide/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_meso_GP',(/vid/),(/gp_mkemeso/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_mean_GP',(/vid/),(/gp_mkemean/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_full_CM',(/vid/),(/cm_mkefull/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_tot_CM' ,(/vid/),(/cm_mketot /),start=(/vertind/))
   call ncinqputvar(ncid,'mke_tide_CM',(/vid/),(/cm_mketide/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_meso_CM',(/vid/),(/cm_mkemeso/),start=(/vertind/))
   call ncinqputvar(ncid,'mke_mean_CM',(/vid/),(/cm_mkemean/),start=(/vertind/))

   call ncputatt(ncid,'mke_full_GP','comment', &
      'Unfiltered RMS of mean kinetic energy GP')
   call ncputatt(ncid,'mke_full_CM','comment', &
      'Unfiltered RMS of mean kinetic energy CM')
   call ncputatt(ncid,'mke_tot_GP','comment', &
      'Total RMS of mean kinetic energy GP')
   call ncputatt(ncid,'mke_tot_CM','comment', &
      'Total RMS of mean kinetic energy CM')
   call ncputatt(ncid,'mke_meso_GP','comment',&
      'Mesoscale RMS of mean kinetic energy GP')
   call ncputatt(ncid,'mke_meso_CM','comment', &
      'Mesoscale RMS of mean kinetic energy CM')
   call ncputatt(ncid,'mke_tide_GP','comment','Tidal RMS of mean kinetic energy GP')
   call ncputatt(ncid,'mke_tide_CM','comment','Tidal RMS of mean kinetic energy CM')
   call ncputatt(ncid,'mke_mean_GP','comment','RMS mke of mean uv GP')
   call ncputatt(ncid,'mke_mean_CM','comment','RMS mke of mean uv CM')
   call ncerr(nf90_close(ncid))

end subroutine mkestat
end module m_mkestat
