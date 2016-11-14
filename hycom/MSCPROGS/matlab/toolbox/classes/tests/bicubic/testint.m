

g=[0 0 8 5];
gy=[0 0 0 0];
gx=[.1 .1 .1 .1];
gxy=[0 0 0 0];

u=0:.10:1.10
t=0:.05:1.05

% Bilinear interpolation
solbl=bilin(g,t,u);


% Bicubic interpolation
c=bicubiccoeffs(g,gx,gy,gxy);
solbc=bicubic(c,t,u);

subplot(1,3,1) ; contourf(solbl',10) ; colorbar
subplot(1,3,2) ; contourf(solbc',10) ; colorbar
subplot(1,3,3) ; contourf(solbc'-solbl',10) ; colorbar
%subplot(1,3,1) ; pcolor(solbl) ; colorbar
%subplot(1,3,2) ; pcolor(solbc) ; colorbar
%subplot(1,3,3) ; pcolor(solbc-solbl) ; colorbar

