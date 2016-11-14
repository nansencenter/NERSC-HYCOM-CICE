  subroutine sort(n,ra,indx)
  implicit none

  integer,  intent(in) ::n
  real,     dimension(n) , intent(inout) :: ra
  integer,  dimension(n),  intent(inout) :: indx

  integer :: l,ir,i,j,indx2
  real :: rra
   
   l=n/2+1
   ir=n

10 continue
   if(l.gt.1)then
     l=l-1
     rra=ra(l)
     indx2=indx(l)
  else
     rra=ra(ir)
     indx2=indx(ir)
     ra(ir)=ra(1)
     indx(ir)=indx(1)
     ir=ir-1
     if(ir.eq.1)then
       ra(1)=rra
       indx(1)=indx2
       return
     endif
   endif

    i=l
    j=l+l
20  if(j.le.ir)then
      if(j.lt.ir)then
        if(ra(j).lt.ra(j+1))j=j+1
      endif
      if(rra.lt.ra(j))then
        ra(i)=ra(j)
        indx(i)=indx(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
    go to 20
    endif

    ra(i)=rra
    indx(i)=indx2
    go to 10
  end subroutine sort
