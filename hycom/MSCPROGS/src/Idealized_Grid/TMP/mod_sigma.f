

      module mod_sigma
      real, parameter :: thref =   1.0e-3
      real, parameter ::       qthref = 1.0/thref

      contains

      real function sig0(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include '../stmt_fns_SIGMA0.h'
      sig0=sig(tt,ss)
      end function


      real function sig2(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include '../stmt_fns_SIGMA2.h'
      sig2=sig(tt,ss)
      end function

      real function sig4(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include '../stmt_fns_SIGMA4.h'
      sig4=sig(tt,ss)
      end function

C-----

      real function kappaf0(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include '../stmt_fns_SIGMA0.h'
      kappaf0=kappaf(tt,ss,pp)
      end function

      real function kappaf2(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include '../stmt_fns_SIGMA2.h'
      kappaf2=kappaf(tt,ss,pp)
      end function

      real function kappaf4(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include '../stmt_fns_SIGMA4.h'
      kappaf4=kappaf(tt,ss,pp)
      end function

C --------

      real function tofsig0(rr,ss)
      implicit none
      real, intent(in) :: rr,ss
      include '../stmt_fns_SIGMA0.h'
      tofsig0=tofsig(rr,ss)
      end function


      real function tofsig2(rr,ss)
      implicit none
      real, intent(in) :: rr,ss
      include '../stmt_fns_SIGMA2.h'
      tofsig2=tofsig(rr,ss)
      end function

      real function tofsig4(rr,ss)
      implicit none
      real, intent(in) :: rr,ss
      include '../stmt_fns_SIGMA4.h'
      tofsig4=tofsig(rr,ss)
      end function

      end module mod_sigma
