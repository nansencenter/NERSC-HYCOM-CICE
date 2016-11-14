function varargout = NCget(filename,varargin)
%% CALL: [v1,v2,...] = NCget(filename,x1,x2,...)
%% where filename is a string
%% and the x1, x2,... have form x={str,rank}
%% where str is a string with the variable name
%% and rank is the rank of the matrix for that variable;

nc          = netcdf(filename);
var_count   = length(varargin);
for j=1:var_count
   x     = varargin{j};
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
   varargout(j)   = {y};
end

close(nc);
