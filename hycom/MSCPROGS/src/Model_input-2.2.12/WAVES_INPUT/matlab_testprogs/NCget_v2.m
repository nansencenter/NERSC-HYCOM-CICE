function [y,y_sf,y_off] = NCget_v2(filename,x)
%% CALL: [v1,v2,...] = NCget(filename,x1,x2,...)
%% where filename is a string
%% and x={str,rank}
%% where str is a string with the variable name
%% and rank is the rank of the matrix for that variable;

nc = netcdf(filename);

str   = x{1};%string with name of variable
rnk   = x{2};%rank of variable
if rnk==1
   y  = nc{ str,1 }(:);
elseif rnk==2
   y  = nc{ str,1 }(:,:);
elseif rnk==3
   y  = nc{ str,1 }(:,:,:);
elseif rnk==4
   y  = nc{ str,1 }(:,:,:,:);
elseif rnk==5
   y  = nc{ str,1 }(:,:,:,:,:);
end

y_sf        = nc{ str,1 }.scale_factor(:);
y_off       = nc{ str,1 }.add_offset(:);
%y_missing   = nc{ str,1 }._FillValue(:);

close(nc);
