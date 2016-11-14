function nc_padheader ( ncfile, num_bytes )
% NC_PADHEADER:  pads the metadata header of a netcdf file
%
% When a netCDF file gets very large, adding new attributes can become
% a time-consuming process.  This can be mitigated by padding the 
% netCDF header with additional bytes.  Subsequent new attributes will
% not result in long time delays unless the length of the new 
% attribute exceeds that of the header.
%
% USAGE:  nc_padheader ( ncfile, num_bytes );
%
% In case of an error, an exception is thrown.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_padheader.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,2,nargin);

[ncid,status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_PADHEADER:MEXNC:OPEN', ncerr );
end

status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_PADHEADER:MEXNC:REDEF', ncerr );
end

%
% Sets the padding to be "num_bytes" at the end of the header section.  
% The other values are default values used by "ENDDEF".
status = mexnc ( '_enddef', ncid, num_bytes, 4, 0, 4 );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_PADHEADER:MEXNC:_ENDDEF', ncerr );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_PADHEADER:MEXNC:CLOSE', ncerr );
end
