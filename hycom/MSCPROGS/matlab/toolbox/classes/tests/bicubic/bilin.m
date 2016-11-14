function sol=bilin(val,t,u)
nt=prod(size(t));
nu=prod(size(u));
for i=1:nt
for j=1:nu
   sol(i,j)=f(val,t(i),u(j));
end
end


% Return coeffs for f(x,y)
function f=f(val,t,u)
   f=(1-t)*(1-u)*val(1) + ...
     (1-t)*(  u)*val(3) + ...
     (  t)*(1-u)*val(2) + ...
     (  t)*(  u)*val(4) ;
