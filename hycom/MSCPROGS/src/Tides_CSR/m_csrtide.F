      module m_csrtide
      contains
      subroutine csrtide (dmjd,dlat,dlon,tide,tlp,tds,isdata,
     .                    u, v, pseudo, radiation)
      use m_lpeqmt
      use m_tptide

c...To compute diurnal and semidiurnal ocean tides from an "orthoweight"
c      file and long period tides from equilibrium theory for a given
c      time, latitude and longitude.

c...Note: An "orthoweight" file can contain geocentric tide, ocean tide,
c         load tide, or tide model difference orthoweights.  If the 
c         orthoweight file in use contains differences or loads then
c         the user should beware of using the "tide" return variable
c         and should use "tds" instead.

c...Coded by Richard Eanes, University of Texas Center for Space Research
c            eanes@csr.utexas.edu

c...Version 2.0 - 03 Nov 94
c               - 13 Nov 94, added better documentation on djmd.
c   Version 3.0 - March 1995, new version of tptide allows compact integer
c                 format for orthoweight file

c...Inputs: (double precision)
c   dmjd    Modified Julian Date in Days, number of days (including fraction)
c           since 17 Nov 1858 00:00.  The standard astronomical julian date
c           (jd) is related to dmjd by jd = 2400000.5 + dmjd.  The modified
c           julian date of 1 Jan 1958 00:00 is 36204.  The modified julian
c           date of 1 Jan 1994 12:00 is 49353.5.
c   dlat    North Latitude (in degrees)
c   dlon    East Longitude (in degrees)
c
c...Outputs (double precision)
c   tds     Diurnal and semidiurnal tide value in millimeters
c   tlp     Long period tide value in millimeters (equilibrium theory)
c   tide    Total tide value in millimeters
c   isdata  Logical, true if location is within the domain of the
c                   tide model, false otherwise.  When false the
c                   tide values are returned as zero.

c The following outputs can be used for constructing maps of the tidal 
c     fields at the major diurnal and semidiurnal constituents using
c     subroutine admit2.
c
c   u,v     Diurnal and semidiurnal orthoweights for input location
c           (can be used in subroutine admit to compute harmonic information)
c   pseudo  Logical, switch passed to subroutine admit determining what
c                    variety of orthoweights are on the tide file
c   radiation  Logical, true if orthoweights modified for radiation potential

cGE      implicit double precision (a-h,o-z)

      logical isdata,pseudo,radiation
      dimension u(6),v(6)

      smjd = dmjd*86400.d0
      call lpeqmt (smjd,dlat,tlp)
      call tptide (dlat,dlon,smjd,tds,isdata,u,v,pseudo,radiation)
      if (isdata) then
        tlp = tlp*10.d0        ! convert from cm to mm
        tds = tds*10.d0        ! convert from cm to mm
        tide = tlp + tds
      else
        tlp = 0.d0
        tds = 0.d0
        tide = 0.d0
      endif

      end subroutine
      end module
