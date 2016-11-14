module m_mersea_prepare


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Handle special situation for generation of MERSEA filenames, based on bulletin
! date, forecast, analysis and hindcast bulletins.


   subroutine mersea_prepare(rungen,rtb,rti,rtd,lprocess,lanalysis)
      use mod_year_info
      implicit none

      type(year_info), intent(inout) :: rti,rtd
      type(year_info), intent(  out) :: rtb
      logical,         intent(out  ) :: lprocess, lanalysis
      character(len=3),intent(in)    :: rungen

      real    :: time1, time2, time
      integer :: year1, year2, refyear
      type(year_info) :: rttmp1,rttmp2
      


      print *,'Processing MERSEA file'

     
      if (rungen=="FOR") then

     
         ! This is the Forecast run.. 
         ! Bulletin date is  Initialization time + 7 days. 
         time=rti%idd+7
         refyear=rti%iyy
         call year_day(time,refyear,rtb,'ecmwf')

     
         ! Get actual days/years of dump date and forecast+6 days
         !time1 = rtb%idd+6
         time1 = rtb%idd
         year1 = rtb%iyy
         time2 = rtd%idd
         year2 = rtd%iyy
         call year_day(time1,year1,rttmp1,'ecmwf')
         call year_day(time2,year2,rttmp2,'ecmwf')

     
         !print *,'rt init    :',rti
         !print *,'rt bulletin:',rttmp1
         !print *,'rt forecast:',rttmp2
         ! We only need the 14 day forecast
         !if ( rttmp1%iyy==rttmp2%iyy .and. rttmp1%idd==rttmp2%idd) then
         if ( (rttmp1%iyy==rttmp2%iyy .and. rttmp1%idd<rttmp2%idd)  .or. &
             (rttmp1%iyy<rttmp2%iyy)) then
            lprocess = .true.
            print *,'This is a part of the "Forecast" time series'
         else
            lprocess=.false.
            print *,'This is not a part of the "Forecast" time series -- skipped'
         end if

     

     
      !TOPAZ2 or earlier! else if (rungen=="RTH") then
      else if (rungen=="ENS") then

     
         ! This is the Ensemble Run. Bulletin date is
         ! Initialization time + 7 days. (These days
         ! are the days leading up to next bulletin).
         time=rti%idd+7
         refyear=rti%iyy
         call year_day(time,refyear,rtb,'ecmwf')
         print *,'rt initial :',rti
         print *,'rt bulletin:',rtb
         print *,'rt forecast:',rtd

     
         ! We only need to process data up to the bulletin date
         if (rtb%iyy==rtd%iyy .and. rtd%idd>=rtb%idd) then
            lprocess=.false.
         elseif (rtb%iyy<rtd%iyy ) then
            lprocess=.false.
            print *,'This is not a part of the "To" time series -- skipped'
         else
            lprocess=.true.
            print *,'This is a part of the "To" time series'
         end if

     
         ! A special case... The bulletin-7 dump is the "Analysis"
         ! of "bulletin - 1". This means that rtd%idd=rti%idd
         if (rti%idd==rtd%idd) then
            lanalysis=.true.
         else
            lanalysis=.false.
         end if

     
      end if

     
      print '(a,a4,1x,a2,1x,a2)','Initial (model start) Year, month, day:', &
         rti%cyy, rti%cmm, rti%cdm
      print '(a,a4,1x,a2,1x,a2)','Bulletin (T0)         Year, month, day:', &
         rtb%cyy, rtb%cmm, rtb%cdm
      print '(a,a4,1x,a2,1x,a2)','Forecast (field dump) Year, month, day:', &
         rtd%cyy, rtd%cmm, rtd%cdm

     
   end subroutine mersea_prepare
end module m_mersea_prepare

