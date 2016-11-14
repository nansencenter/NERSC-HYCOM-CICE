function bool = nc_iscoordvar ( ncfile, varname )
% NC_ISCOORDVAR:  yes if the given variable is also a coordinate variable.
%
% A coordinate variable is a variable with just one dimension.  That 
% dimension has the same name as the variable itself.
%
% USAGE:  bool = nc_iscoordvar ( ncfile, varname );
%
% PARAMETERS:
% Input:
%     ncfile:  
%        Input netcdf file name.
%     varname:  
%        variable to check
% Output:
%     bool:
%         1 if the variable is a coordinate variable
%         0 if the variable is not a coordinate variable
%
% Throws an exception in case of an error condition.
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_iscoordvar.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


snc_nargchk(2,2,nargin);
snc_nargoutchk(1,1,nargout);

%
% Assume that the answer is no until we know that it is yes.

bool = false;

ncvar = nc_getvarinfo ( ncfile, varname );

%
% Check that it's not a singleton.  If it is, then the answer is no.
if isempty(ncvar.Dimension)
	bool = false;
	return
end

%
% Check that the names are the same.
if strcmp ( ncvar.Dimension{1}, varname )
	bool = true;
end


return;

