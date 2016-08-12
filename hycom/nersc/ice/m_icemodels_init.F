#if defined (ICESTATE) && defined (ICE)
#error CPP FLAGS You have defined both ICE and ICESTATE
#endif
      module m_icemodels_init
C --- This routine takes care of all ice model initialization 
C --- giving a CLEANER HYCOM code.
      contains
      subroutine icemodels_init()
      use mod_xc
#if defined (ICESTATE)
      use mod_icestate_init
#endif
      use mod_common_ice
      implicit none
      include 'common_blocks.h'

c
c --- Initialize variables in mod_common_ice
      call icedat
c
c --- Initialize Hunke&Dukowicz EVP dynamics model
#if defined (EVP)
      if (mnproc==1) write(lp,'(a)') 'Initializing EVP from hycom'
      call evp_init_from_hycom()
      if (mnproc==1) write(lp,'(a)') 'Initializing EVP from input file'
      call evp_init()
#endif
c
c --- Initialize ICESTATE thermodynamic/kinematic model
#if defined (ICESTATE)
      if (mnproc==1) write(lp,'(a)') 'ICESTATE init:'
      call icestate_read_infile
      call icestate_init(baclin,depths,plat*radian,thref)
#endif
      end subroutine icemodels_init
      end module m_icemodels_init
