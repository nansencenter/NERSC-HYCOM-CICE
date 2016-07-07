          module globals
          implicit none

	  REAL*4, ALLOCATABLE :: A(:,:,:),p(:,:,:),AZ(:,:,:),wk1(:,:,:)
	  REAL*4, ALLOCATABLE :: lon(:,:),lat(:,:),A2d(:,:)
          REAL*4, ALLOCATABLE :: depth(:,:),ubaro(:,:),vbaro(:,:)
	  REAL*4, ALLOCATABLE :: dp(:,:,:)
          INTEGER, DIMENSION(20) ::layers
	  REAL*4, ALLOCATABLE      :: sigma(:)
          REAL*4                    :: filler,SPVAL
	  LOGICAL                 :: layer,lyr2z
          CHARACTER(LEN=5)        :: Vaxis
          CHARACTER(LEN=10)       :: lon_step,lat_step
          INTEGER   :: Time,IDM,JDM,kdm,oldmode,field_julian_date
          CHARACTER(LEN=100) :: name,
     &                          long_name,
     &                          units,
     &                          fname,
     &                          Title,
     &                          Institution,
     &                          history,
     &                          EXPT,
     &                          ARTYPE,
     &                          easting,
     &                          northing,
     &                          reference,
     &                          domain,
     &                          domain_name,
     &                          field_type,
     &                          field_date,
     &                          julian_date_unit,
     &                          forecast_range,
     &                          forecast_type,
     &                          bulletin_type,
     &                          bulletin_date,
     &                          standard_name,
     &                          unit_long,
     &                          comments,
     &                          conventions,
     &                          source,
     &                          lon_standard_name,
     &                          lat_standard_name,
     &                          depth_standard_name,
     &                          lon_axis,
     &                          lat_axis,
     &                          depth_axis,
     &                          z_positive_direction,
     &                          data_type,
     &                          longname
 
          end module globals
