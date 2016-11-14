function gpdens_exceedence(infile)
%function gp_timeseries(infile)
%
% function plots filtered, unfiltered and tidal timeseries

vlev = nc_varget ( infile, 'vlevel' )
nvlev = prod(size(vlev));
elev = nc_varget ( infile, 'e_levels' )
nelev = prod(size(elev));


e_GP= nc_varget ( infile, 'exceedence_GP' );

whos
cvals=[ 'r' ; 'g' ; 'b' ; 'c' ; 'm' ; 'k' ; 'y' ; 'g'];

[tmp,I] = sort(vlev);
figure(1); clf; hold on;
for i=1:nelev
   plot(e_GP(I,i),-vlev(I),'LineWidth',2,'Color',cvals(i))
   %plot(e_GP(I,i),-vlev(I),'LineWidth',2)
end
xlabel('cm/s')
ylabel('Depth')

