module mod_types
   type fieldheaders
      character(len=8) fieldname
      integer          layer
      integer          ix
      integer          iy
      integer          ia
      integer          ib
      integer          ic
      integer          length
   end type fieldheaders

   type fields
      character(len=8) fextract
      integer          laya
      integer          layb
      logical          option
      logical          vecflag
      character(len=8) vecpost
   end type fields
end module mod_types

