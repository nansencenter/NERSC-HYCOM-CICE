        subroutine wrtNcdf(ii,jj,kk,ix1,jy1,fill_value,dim3)

        use TYPESIZES ! NetCDF needs this
        use NETCDF    ! NetCDF fortran 90 interface
        use globals   ! global info arrays

        implicit none

        integer          :: ii,jj,kk,iidm,jjdm,ix1,jy1
         REAL             :: fill_value,missing_value
        integer          :: ncfileId, status, varId
        integer          :: pLatDimId,pLonDimId,pLatVarId,pLonVarId
        integer          :: ptDimId, ptVarId,lyrDimID,lyrVarID
        character        :: ncfile*100,file*100
        LOGICAL          :: dim3,ex,transport
        real             :: hmin(999),hmax(999)


        ncfile= trim(fname) // trim(".nc")
        iidm=ix1+ii-1
        jjdm=jy1+jj-1
        
        missing_value = fill_value 
        
    
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
          call ncheck(nf90_def_dim(ncfileId,"layer",kk,lyrDimId))
          END IF
            

!!!!!!!!!!!!!!!!!!!!!GLOBAL ATTRIBUTES !!!!!!!!!!!!!!!!!!!!!!!!!!!


             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "conventions",
     &                             "CF-1.0"))
             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "title",
     &                              Title))
             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                             "institution",
     &                              Institution))
             
            
            call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "experiment",
     &                               EXPT))
            call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "history",history))
     &                        

        

             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "domain_name",
     &                                domain_name))


             call ncheck(nf90_put_att(ncfileID,nf90_global,
     &                               "field_date",
     &                                field_date))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! create the variables and attributes

           call ncheck(nf90_def_var(ncfileId,"Longitude", nf90_float,
     &      pLonDimId,pLonVarId))
           call ncheck(nf90_put_att(ncfileId,pLonVarId,
     &      "units","degrees_east"))

           call ncrange(lon(ix1:iidm,jy1:jjdm),ii,jj,1, 
     &     fill_value, hmin,hmax)

           call ncheck(nf90_put_att(ncfileID,pLonvarID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))
           
!           call ncheck(nf90_put_att(ncfileID,pLonVarID,
!     &                             "step",
!     &                              lon_step))

           call ncheck(nf90_put_att(ncfileID,pLonVarID,
     &                             "standard_name",
     &                              lon_standard_name))


           call ncheck(nf90_put_att(ncfileID,pLonVarID,
     &                              "axis",lon_axis))
     &            
          call ncheck(nf90_def_var(ncfileId,"Latitude", nf90_float,
     &      pLatDimId,pLatVarId))
          call ncheck(nf90_put_att(ncfileId,pLatVarId,
     &      "units","degrees_north"))


           call ncrange(lat(ix1:iidm,jy1:jjdm),ii,jj,1,
     &     fill_value, hmin,hmax)

           call ncheck(nf90_put_att(ncfileID,pLatvarID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))

!           call ncheck(nf90_put_att(ncfileID,pLatVarID,
!     &                             "step",
!     &                              lat_step))


           call ncheck(nf90_put_att(ncfileID,pLatVarID,
     &                             "standard_name",
     &                              lat_standard_name))


           call ncheck(nf90_put_att(ncfileID,pLatVarID,
     &                             "axis",lat_axis))
     &


          IF(dim3) THEN
          call ncheck(nf90_def_var(ncfileId,"layer", nf90_float,
     &      lyrDimId,lyrVarId))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         END IF ! 3d loop





          

          ! leave define mode
          call ncheck(nf90_enddef(ncfileId))
          

          ! write data into coordinate variables
          call ncheck(nf90_put_var(ncfileID,pLatVarId,lat(1,jy1:jjdm)))
          call ncheck(nf90_put_var(ncfileID,pLonVarId,lon(ix1:iidm,1)))


          IF(dim3) THEN

              call ncheck(nf90_put_var(ncfileID,LyrVarId,sigma(:)))

            END IF
          
           
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
            call ncheck(nf90_inq_dimid(ncfileID,"layer",lyrDimID))
            call ncheck(nf90_def_var(ncfileId,trim(name),nf90_float,
     &        (/pLonDimId, pLatDimId,LyrDimId/), varID))
            call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
     &        fill_value))
            call ncheck(nf90_put_att(ncfileId,varId,"units",units))
            call ncheck(nf90_put_att(ncfileId,varId,"standard_name",
     &        standard_name))

            call ncheck(nf90_put_att(ncfileId,varId,"long_name",
     &        long_name))

            call ncheck(nf90_put_att(ncfileId,varId,"unit_long",
     &        unit_long))

            call ncheck(nf90_put_att(ncfileId,varId,"missing_value",
     &        missing_value))

            call ncheck(nf90_put_att(ncfileID,varID,
     &        "valid_range", (/hmin(1), hmax(1)/)))



          ELSE

            call ncrange(A2D,ii,jj,1, fill_value, hmin,hmax)
            call ncheck(nf90_inq_dimid(ncfileID,"Latitude",pLatDimID))
            call ncheck(nf90_inq_dimid(ncfileID,"Longitude",pLonDimID))
            call ncheck(nf90_def_var(ncfileId,trim(name),nf90_float,
     &       (/pLonDimId, pLatDimId/),varID))
            call ncheck(nf90_put_att(ncfileId,varId,"_FillValue",
     &       fill_value))
            call ncheck(nf90_put_att(ncfileId,varId,"units",
     &       units))
            call ncheck(nf90_put_att(ncfileId,varId,"standard_name",
     &       standard_name))
            call ncheck(nf90_put_att(ncfileId,varId,"long_name",
     &        long_Name))
            call ncheck(nf90_put_att(ncfileId,varId,"unit_long",
     &        unit_long))
            call ncheck(nf90_put_att(ncfileId,varId,"missing_value",
     &        missing_value))
            call ncheck(nf90_put_att(ncfileID,varID,
     &        "valid_range", (/hmin(1), hmax(1)/)))

          END IF



      ! leave define mode
             call ncheck(nf90_enddef(ncfileId))



      IF(dim3) THEN
             call ncheck(nf90_inq_varid(ncfileId,name,varId))




        ! print *, "writing 3d" 

       call ncheck(nf90_put_var(ncfileId,varId,
     &    A(ix1:iidm,jy1:jjdm,1:kk),
     & start=(/1,1,1/),count=(/ii,jj,kk/)))


       
      ELSE

      !print *, "writitng 2d"
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
  
      
