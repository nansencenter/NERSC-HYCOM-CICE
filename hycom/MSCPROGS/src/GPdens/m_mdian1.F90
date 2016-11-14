module m_mdian1
contains
subroutine mdian1(x,n,xmed)
  !use m_sort

  dimension x(n)
  integer sortindx(n)
  call sort(n,x,sortindx)
  n2=n/2

  if(2*n2.eq.n)then
   xmed=0.5*(x(n2)+x(n2+1))
   else
   xmed=x(n2+1)
  endif
  end subroutine mdian1
end module m_mdian1
