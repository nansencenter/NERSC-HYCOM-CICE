        subroutine ncfinit(ii,jj,kk,ix1,jy1,fill_value,dim3,transport)

        use TYPESIZES ! NetCDF needs this
        use NETCDF    ! NetCDF fortran 90 interface
        use globals   ! global info arrays

        implicit none

        integer          :: ii,jj,kk,iidm,jjdm,ix1,jy1
        REAL             :: fill_value
        integer          :: ncfileId, status, varId
        integer          :: pLatDimId,pLonDimId,pLatVarId,pLonVarId
        integer          :: ptDimId, ptVarId,lyrDimID,lyrVarID
        character        :: ncfile*100,file*100
        LOGICAL          :: dim3,ex,transport
        real             :: hmin(999),hmax(999)
        ncfile= trim(fname) // trim(".nc")
        iidm=ix1+ii-1
        jjdm=jy1+jj-1
     
        
    
        INQUIRE (file = ncfile, exist = ex )

          IF ( .not. ex ) then
!       If file does not exist  create a new NetCDF and write data to it

          call ncheck(nf90_create(trim(ncfile),nf90_noclobber,ncfileId))
         

          ! define the dimensions
          ! set nofill mode

          call ncheck(nf90_set_fill(ncfileID,nf90_nofill,oldmode))
          call ncheck(nf90_def_dim(ncfileId,"Longitude",ii,pLonDimId))
          call ncheck(nf90_def_dim(ncfileId,"Latitude",jj,pLatDimId))
!         if 3D variable then create additional variable
          IF (dim3) THEN
          call ncheck(nf90_def_dim(ncfileId,vaxis,kk,lyrDimId))
          END IF
            
!          call ncheck(nf90_def_dim(ncfileId,"Time",1,ptDimId))

!!!!!!!!!!!!!!!!!!!!!GLOBAL ATTRIBUTES !!!!!!!!!!!!!!!!!!!!!!!!!!!


             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "Conventions",
     &                             "CF-1.0"))
             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "title",
     &                              Title))
             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "institution",
     &                              Institution))
             
             If(ARTYPE.eq.'daily') THEN          

             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                              "source",
     &                              "HYCOM archive file"))
             else
   
             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                              "source",
     &                              "HYCOM mean archive file"))
             endif
            
            call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "experiment",
     &                               EXPT))
            call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "history","a2data"))
     &                        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




          ! create the variables and attributes

          call ncheck(nf90_def_var(ncfileId,"Longitude", nf90_float,
     &      pLonDimId,pLonVarId))
          call ncheck(nf90_put_att(ncfileId,pLonVarId,
     &      "units","degrees_east"))
          call ncheck(nf90_def_var(ncfileId,"Latitude", nf90_float,
     &      pLatDimId,pLatVarId))
          call ncheck(nf90_put_att(ncfileId,pLatVarId,
     &      "units","degrees_north"))

          IF(dim3) THEN
          call ncheck(nf90_def_var(ncfileId,vaxis, nf90_int,
     &      lyrDimId,lyrVarId))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (vaxis .eq. 'depth') THEN

          call ncheck(nf90_put_att(ncfileId,lyrVarId,
     &      "units","m"))
          ELSE
          call ncheck(nf90_put_att(ncfileId,LyrVarId,
     &      "units","layer"))
          print *, "layer"
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          END IF





!          call ncheck(nf90_def_var(ncfileId,"Time", nf90_float,
!     &      ptDimId,ptVarId))
!          call ncheck(nf90_put_att(ncfileId,ptVarId,
!     &      "units","days since 1900-12-31"))

!          call ncheck(nf90_put_att(ncfileId,ptVarId,
!     &      "time_origin","1900-12-31"))
          

          ! leave define mode
          call ncheck(nf90_enddef(ncfileId))
          

          ! write data into coordinate variables
          call ncheck(nf90_put_var(ncfileID,pLatVarId,lat(1,jy1:jjdm)))
          call ncheck(nf90_put_var(ncfileID,pLonVarId,lon(ix1:iidm,1)))


          IF(dim3) THEN

            IF (vaxis .eq. 'depth') THEN
              call ncheck(nf90_put_var(ncfileID,LyrVarId,depths(1:kk)))
            ELSE
              call ncheck(nf90_put_var(ncfileID,LyrVarId,layers(1:kk)))
            END IF

            END IF
          
!            call ncheck(nf90_put_var(ncfileId,ptvarId,time))
           
            call ncheck(nf90_redef(ncfileID))
          

          ELSE
          

          ! model variables
c*****************************************************************          
          ! Write data to existing NC file

          call ncheck(nf90_open(trim(ncfile),nf90_write,ncfileID))
          call ncheck(nf90_redef(ncfileID))

          END IF
                
       
          IF (dim3) THEN
            call ncrange(A,ii,jj,kk, fill_value, hmin,hmax)   
            call ncheck(nf90_inq_dimid(ncfileID,"Latitude",pLatDimID))
            call ncheck(nf90_inq_dimid(ncfileID,"Longitude",pLonDimID))
            call ncheck(nf90_inq_dimid(ncfileID,vaxis,lyrDimID))
!            call ncheck(nf90_inq_dimid(ncfileID,"Time",ptDimId))
    
          call ncheck(nf90_def_var(ncfileId,name,nf90_float,
     &        (/pLonDimId, pLatDimId,lyrDimId/), varID))
          call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
     &        fill_value))
          call ncheck(nf90_put_att(ncfileId,varId,"units",units))
          call ncheck(nf90_put_att(ncfileId,varId,"standard_name",
     &        longName))
          call ncheck(nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))

!!!!!!!!!!!!!!!!!!!!!!!! Transport!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!          IF(name.eq.'uvel' .and. transport) THEN
!          call ncheck(nf90_def_var(ncfileId,"udp",nf90_float,
!     &          (/pLonDimId, pLatDimId,ptDimId/),varID))
!          call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
!     &           fill_value))
!
!          call ncheck(nf90_put_att(ncfileId,varId,"units",
!     &      "m2s-1"))
!          call ncheck(nf90_put_att(ncfileId,varId,"long_name",
!     &      "Depth_Integrated_Eastward_Sea_Water_Velocity"))

!          ELSE IF(name.eq.'vvel' .and. transport) THEN

!            call ncheck(nf90_def_var(ncfileId,"vdp",nf90_float,
!     &       (/pLonDimId, pLatDimId,ptDimId/),varID))
!            call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
!     &       fill_value))

!            call ncheck(nf90_put_att(ncfileId,varId,"units",
!     &        "m2s-1"))
!            call ncheck(nf90_put_att(ncfileId,varId,"long_name",
!     &        "Depth_Integrated_Northward_Sea_Water_Velocity"))

!            END IF
          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!Transport!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          ELSE
        call ncrange(A2D,ii,jj,1, fill_value, hmin,hmax)
        call ncheck(nf90_inq_dimid(ncfileID,"Latitude",pLatDimID))
        call ncheck(nf90_inq_dimid(ncfileID,"Longitude",pLonDimID))
!        call ncheck(nf90_inq_dimid(ncfileID,"Time",ptDimId))
          
       

          call ncheck(nf90_def_var(ncfileId,trim(name),nf90_float,
     &       (/pLonDimId, pLatDimId/),varID))
          call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
     &       fill_value))
         
          call ncheck(nf90_put_att(ncfileId,varId,"units",
     &       units))
          call ncheck(nf90_put_att(ncfileId,varId,"standard_name",
     &       longName))
          
          call ncheck(nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))

          END IF



      ! leave define mode
      call ncheck(nf90_enddef(ncfileId))



      IF(dim3) THEN
!      call ncheck(nf90_inq_varid(ncfileId,name,varId))
!      IF (vaxis .eq. 'depth') THEN

!!!!!!!!!!!!!!!!!!!! If Transport!!!!!!!!!!!!!!!!!!!!!
!      IF(name.eq.'uvel' .and. transport) THEN

        
!      call ncheck(nf90_put_var(ncfileId,varId,
!     & AZ(ix1:iidm,jy1:jjdm,1:kk),
!     & start=(/1,1,1,1/),count=(/ii,jj,kk,1/)))
!      call ncheck(nf90_inq_varid(ncfileId,"udp",varId))
!      call ncheck(nf90_put_var(ncfileId,varId,
!     &    A2d(ix1:iidm,jy1:jjdm),
!     &  start=(/1,1,1/),count=(/ii,jj,1/)))

      
!      ELSE IF(name .eq. 'vvel' .and. transport) THEN

!       call ncheck(nf90_put_var(ncfileId,varId,
!     &   AZ(ix1:iidm,jy1:jjdm,1:kk),
!     &   start=(/1,1,1,1/),count=(/ii,jj,kk,1/)))
!       call ncheck(nf90_inq_varid(ncfileId,"vdp",varId))
!       call ncheck(nf90_put_var(ncfileId,varId,
!     &   A2d(ix1:iidm,jy1:jjdm),
!     &   start=(/1,1,1/),count=(/ii,jj,1/)))
!      END IF

!!!!!!!!!!!!!!!!!!!!!! End Transport!!!!!!!!!!!!!!!!!!!!!!

!      ELSE
!!! Write layer
         print *, "writing 3d" 

       call ncheck(nf90_put_var(ncfileId,varId,
     &    A(ix1:iidm,jy1:jjdm,1:kk),
     & start=(/1,1,1/),count=(/ii,jj,kk/)))


!      ENDIF
       
      ELSE

!! Write 2d
      print *, "writitng 2d"
      call ncheck(nf90_inq_varid(ncfileId,name,varId))
      call ncheck(nf90_put_var(ncfileId,varId,
     &  A2d(ix1:iidm,jy1:jjdm),
     &  start=(/1,1/),count=(/ii,jj/)))
      END IF



      ! close file
      call ncheck(nf90_close(ncfileId))


      END






        subroutine ncheck(status)
        use netcdf   ! NetCDF fortran 90 interface
        implicit none
        integer, intent(in) :: status

c       subroutine to handle NetCDF errors

       if (status /= nf90_noerr) then
          write(*,*)   'error in horout - from NetCDF library'
          write(*,*)   trim(nf90_strerror(status))
          stop
        end if
        end subroutine ncheck
  
      
