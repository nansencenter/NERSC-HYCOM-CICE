      module mod_sigma
      real, parameter :: thref =   1.0e-3
      real, parameter ::       qthref = 1.0/thref

      contains

      real function sig0(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include 'stmt_fns_SIGMA0.H'
      sig0=sig(tt,ss)
      end function


      real function sig2(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include 'stmt_fns_SIGMA2.H'
      sig2=sig(tt,ss)
      end function

      real function sig4(tt,ss)
      implicit none
      real, intent(in) :: tt,ss
      include 'stmt_fns_SIGMA4.H'
      sig4=sig(tt,ss)
      end function

C-----

      real function kappaf0(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include 'stmt_fns_SIGMA0.H'
      kappaf0=kappaf(tt,ss,pp)
      end function

      real function kappaf2(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include 'stmt_fns_SIGMA2.H'
      kappaf2=kappaf(tt,ss,pp)
      end function

      real function kappaf4(tt,ss,pp)
      implicit none
      real, intent(in) :: tt,ss,pp
      include 'stmt_fns_SIGMA4.H'
      kappaf4=kappaf(tt,ss,pp)
      end function

      end module mod_sigma
