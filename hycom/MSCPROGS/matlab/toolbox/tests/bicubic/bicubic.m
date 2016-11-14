function sol=bicubic(coeffs,t,u)
%f  (x,y) = sum_i=1^4 sum_j=1,4 c_ij t^(i-1) t^(j-1)
%fx (x,y) = sum_i=2^4 sum_j=1,4 c_ij (i-1) * t^(i-2) t^(j-1)
%fy (x,y) = sum_i=1^4 sum_j=2,4 c_ij (j-1) * t^(i-1) t^(j-2)
%fxy(x,y) = sum_i=2^4 sum_j=2,4 c_ij (j-1) * (j-1) * t^(i-2) t^(j-2)

coeffs=reshape(coeffs,4,4);
nt=prod(size(t));
nu=prod(size(u));
sol=zeros(size(nt,nu));

for i=1:nt
for j=1:nu
   sol(i,j)=f(coeffs,t(i),u(j));
end
end


% Return coeffs for f(x,y)
function f=f(c,t,u)
   f=0.;
   for i=1:4
   for j=1:4
      f=f+c(i,j) * t^(i-1) * u^(j-1);
   end 
   end
