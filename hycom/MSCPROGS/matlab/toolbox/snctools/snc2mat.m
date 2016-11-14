function snc2mat ( ncfile, matfile )
% SNC2MAT:  saves netcdf file to *.mat format.  
%
% SNC2MAT(NCFILE,MATFILE) will save the netCDF file NCFILE to the mat-file 
% MATFILE.  This function is deprecated and may disappear in a future release
% of SNCTOOLS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: snc2mat.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wid = sprintf ( 'SNCTOOLS:%s:deprecatedMessage', lower(mfilename) );
msg = sprintf( '%s is deprecated and may be removed in a future version of SNCTOOLS.', upper(mfilename) );
warning ( wid, msg );



%
% create the MATLAB file
ncdata = nc_getall ( ncfile );

fnames = fieldnames ( ncdata );
save_command = '';
global_atts = [];
for j = 1:length(fnames)
	theVar = fnames{j};
	if ( strcmp(theVar,'global_atts' ) )
		global_atts = ncdata.global_atts;
	else
		command = sprintf ( '%s = ncdata.%s;', theVar, theVar );
		eval(command);
		save_command = sprintf ( '%s''%s'',', save_command, theVar );
	end
end
if ~isempty(global_atts)
	save_command = sprintf ( '%s''global_atts''', save_command );
else
	%
	% This chops off a bad comma that's not needed if no global attributes.
	save_command(end) = '';
end
save_command = sprintf ( 'save ( matfile, %s );', save_command );
try
	eval(save_command);
catch
	msg = sprintf ( 'Could not execute ''%s'', got error message ''%s''\n', save_command, lasterr );
	snc_error ( 'SNCTOOLS:snc2mat:badCommand', msg );
end


	

return
