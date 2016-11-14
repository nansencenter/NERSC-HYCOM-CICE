function nc_addhist ( ncfile, attval )
% NC_ADDHIST:  adds text to a global history attribute
%
% NC_ADDHIST(NCFILE,TEXT) adds the TEXT string to the standard convention
% "history" global attribute of the netCDF file NCFILE.  The string is 
% prepended, rather than appended.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_addhist.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,2,nargin);

if ~exist(ncfile,'file')
	snc_error ( 'SNCTOOLS:NC_ADDHIST:badFilename', '%s does not exist', ncfile );
end
if ~ischar(attval)
	snc_error ( 'SNCTOOLS:NC_ADDHIST:badDatatype', ...
	        'history attribute addition must be character.' );
end


try
	old_hist = nc_attget ( ncfile, nc_global, 'history' );
catch
	%
	% The history attribute must not have existed.  That's ok.
	old_hist = '';
end


if isempty(old_hist)
	new_history = sprintf ( '%s:  %s', datestr(now), attval );
else
	new_history = sprintf ( '%s:  %s\n%s', datestr(now), attval, old_hist );
end
nc_attput ( ncfile, nc_global, 'history', new_history );



