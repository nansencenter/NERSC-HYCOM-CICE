function nc_add_dimension ( ncfile, dimension_name, dimension_length )
% NC_ADD_DIMENSION:  adds a dimension to an existing netcdf file
%
% USAGE:  nc_add_dimension ( ncfile, dimension_name, dimension_size );
%
% PARAMETERS:
% Input:
%     ncfile:  path to netcdf file
%     dimension_name:  name of dimension to be added
%     dimension_size:  length of new dimension.  If zero, it will be an
%         unlimited dimension.
% Output:
%     none
%
% In case of an error, an exception is thrown.
%
% Because of underlying limitations, this m-file requires mexnc.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_add_dimension.m 2305 2007-08-31 15:24:15Z johnevans007 $
% $LastChangedDate: 2007-08-31 11:24:15 -0400 (Fri, 31 Aug 2007) $
% $LastChangedRevision: 2305 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(3,3,nargin);

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if status
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_DIMENSION:openFailed';
	snc_error ( error_id, ncerr );
end

status = mexnc ( 'redef', ncid );
if status
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_DIMENSION:redefFailed';
	snc_error ( error_id, ncerr );
end

[dimid, status] = mexnc ( 'def_dim', ncid, dimension_name, dimension_length );
if status
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_DIMENSION:defdimFailed';
	snc_error ( error_id, ncerr );
end

status = mexnc ( 'enddef', ncid );
if status
	mexnc ( 'close', ncid );
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_DIMENSION:enddefFailed';
	snc_error ( error_id, ncerr );
end


status = mexnc ( 'close', ncid );
if status 
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_DIMENSION:closeFailed';
	snc_error ( error_id, ncerr );
end



return












