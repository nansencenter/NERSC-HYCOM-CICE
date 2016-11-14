function create_empty_file ( ncfile )
% CREATE_EMPTY_FILE:  Does just that, makes an empty netcdf file.
%
% USAGE:  create_empty_file ( ncfile );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: create_empty_file.m 1781 2006-09-08 11:12:02Z johnevans007 $
% $LastChangedDate: 2006-09-08 07:12:02 -0400 (Fri, 08 Sep 2006) $
% $LastChangedRevision: 1781 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ok, first create the first file
[ncid_1, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end

%
% CLOSE
status = mexnc ( 'close', ncid_1 );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end
return
