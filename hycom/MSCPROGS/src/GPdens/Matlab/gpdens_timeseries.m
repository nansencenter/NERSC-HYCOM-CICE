function gpdens_timeseries(infile,choice)
%function gp_timeseries(infile)
%
% function plots filtered, unfiltered and tidal timeseries

if (strcmp(choice,'spd'))
   var_gp_U = nc_varget ( infile, 'spd_GP_U' );
   var_gp_F = nc_varget ( infile, 'spd_GP_F' );
   var_gp_T = nc_varget ( infile, 'spd_GP_T' );
elseif (strcmp(choice,'u'))
   var_gp_U = nc_varget ( infile, 'u_GP_U' );
   var_gp_F = nc_varget ( infile, 'u_GP_F' );
   var_gp_T = nc_varget ( infile, 'u_GP_T' );
elseif (strcmp(choice,'v'))
   var_gp_U = nc_varget ( infile, 'v_GP_U' );
   var_gp_F = nc_varget ( infile, 'v_GP_F' );
   var_gp_T = nc_varget ( infile, 'v_GP_T' );
end
vlev = nc_varget ( infile, 'vlevel' );



[tmp,I] = sort(vlev);

figure(1); clf
nlev=prod(size(I));
for i=1:nlev
   subplot(nlev,1,i) ; 
   hold on;
   plot(var_gp_U(I(i),:),'Color','b','LineWidth',2)
   plot(var_gp_F(I(i),:),'Color','r','LineWidth',2)
   plot(var_gp_T(I(i),:),'Color','g','LineWidth',2)
end



