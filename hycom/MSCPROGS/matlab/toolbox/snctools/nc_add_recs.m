function nc_add_recs ( ncfile, new_data, varargin )
% NC_ADD_RECS:  add records onto the end of a netcdf file
%
% USAGE:  nc_add_recs ( ncfile, new_data, unlimited_dimension );
% 
% INPUT:
%   ncfile:  netcdf file
%   new_data:  Matlab structure.  Each field is a data array
%      to be written to the netcdf file.  Each array had
%      better be the same length.  All arrays are written
%      in the same fashion.
%   unlimited_dimension:
%      Optional.  Name of the unlimited dimension along which the data 
%      is written.  If not provided, we query for the first unlimited 
%      dimension (looking ahead to HDF5/NetCDF4).
%     
% OUTPUT:
%   None.  In case of an error, an exception is thrown.
%
% AUTHOR: 
%   johnevans@acm.org
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_add_recs.m 2309 2007-08-31 20:30:56Z johnevans007 $
% $LastChangedDate: 2007-08-31 16:30:56 -0400 (Fri, 31 Aug 2007) $
% $LastChangedRevision: 2309 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


snc_nargchk(2,3,nargin);


%
% Check that we were given good inputs.
if ~isstruct ( new_data )
	err_id = 'SNCTOOLS:NC_ADD_RECS:badStruct';
	snc_error ( err_id, '2nd input argument must be a structure .\n' );
end

%
% Check that each field of the structure has the same length.
varnames = fieldnames ( new_data );
num_fields = length(varnames);
if ( num_fields <= 0 )
	err_id = 'SNCTOOLS:NC_ADD_RECS:badRecord';
	snc_error ( err_id, 'data record cannot be empty' );
end
field_length = zeros(num_fields,1);
for j = 1:num_fields
	command = sprintf ( 'field_length(j) = size(new_data.%s,1);', varnames{j} );
	eval ( command );
end
if any(diff(field_length))
	err_id = 'SNCTOOLS:NC_ADD_RECS:badFieldLengths';
	snc_error ( err_id, 'Some of the fields do not have the same length.\n' );
end

%
% So we have this many records to write.
record_count(1) = field_length(1);


[unlim_dimname, unlim_dimlen, unlim_dimid] = get_unlimdim_info ( ncfile, varargin{:} );

varsize = get_all_varsizes ( ncfile, new_data, unlim_dimid );


%
% So we start writing here.
record_corner(1) = unlim_dimlen;



%
% write out each data field, as well as the minimum and maximum
input_variable = fieldnames ( new_data );
num_vars = length(input_variable);
for i = 1:num_vars

	current_var = input_variable{i};
	%fprintf ( 1, '%s:  processing %s...\n', mfilename, current_var );

	current_var_data = new_data.(current_var);
	var_buffer_size = size(current_var_data);

	netcdf_var_size = varsize.(current_var);

	corner = zeros( 1, length(netcdf_var_size) );
	count = ones( 1, length(netcdf_var_size) );

	corner(1) = record_corner(1);
	count(1) = record_count(1);


	
	for j = 2:numel(var_buffer_size)
		if ( var_buffer_size(j) > 1 )
			count(j) = var_buffer_size(j);
		end
	end



	%
	% Ok, we are finally ready to write some data.
	nc_varput ( ncfile, current_var, current_var_data, corner, count );
   

end


return









function varsize = get_all_varsizes ( ncfile, new_data,unlimited_dimension_dimid )

[ncid,status ]=mexnc( 'open', ncfile, nc_nowrite_mode );
if status ~= 0
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_RECS:openFailed';
	snc_error ( error_id, ncerr );
end



%
% For each field of "new_data" buffer, inquire as to the dimensions in the
% NetCDF file.  We need this data to properly tell nc_varput how to write
% the data
input_variable = fieldnames ( new_data );
num_vars = length(input_variable);
varsize = [];
for j = 1:num_vars

	[varid, status] = mexnc('INQ_VARID', ncid, input_variable{j} );
	if ( status ~= 0 )
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		error_id = 'SNCTOOLS:NC_ADD_RECS:inq_varidFailed';
		snc_error ( error_id, ncerr );
	end

	[dimids, status] = mexnc('INQ_VARDIMID', ncid, varid);
	if ( status ~= 0 )
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		error_id = 'SNCTOOLS:NC_ADD_RECS:inq_vardimidFailed';
		snc_error ( error_id, ncerr );
	end
	ndims = length(dimids);
	dimsize = zeros(ndims,1);


	%
	% make sure that this variable is defined along the unlimited dimension.
	if ~any(find(dimids==unlimited_dimension_dimid))
		mexnc('close',ncid);
		format = 'variable %s must be defined along unlimited dimension %s.\n';
		snc_error ( 'SNCTOOLS:NC_ADD_RECS:missingUnlimitedDimension', ...
		        format, input_variable{j}, unlimited_dimension_name );
	end

	for k = 1:ndims
		[dim_length, status] = mexnc('INQ_DIMLEN', ncid, dimids(k) );
		if ( status ~= 0 )
			mexnc('close',ncid);
			ncerr = mexnc ( 'strerror', status );
			error_id = 'SNCTOOLS:NC_ADD_RECS:inq_dimlenFailed';
			snc_error ( error_id, ncerr );
		end
		dimsize(k) = dim_length;
	end

	varsize.(input_variable{j}) = dimsize;

end

status = mexnc('close',ncid);
if status ~= 0 
	ncerr = mexnc ( 'strerror', status );
	error_id = 'SNCTOOLS:NC_ADD_RECS:closeFailed';
	snc_error ( error_id, ncerr );
end







function [dimname, dimlen, dimid] = get_unlimdim_info ( ncfile, varargin )

[ncid,status ]=mexnc( 'open', ncfile, nc_nowrite_mode );
if status ~= 0
	mexnc('close',ncid);
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_ADD_RECS:openFailed', ncerr );
end



%
% If we were not given the name of an unlimited dimension, get it now
if nargin < 2
	[dimid, status] = mexnc ( 'inq_unlimdim', ncid );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'SNCTOOLS:NC_ADD_RECS:inq_unlimdimFailed', ncerr );
	end

	[dimname, status] = mexnc ( 'INQ_DIMNAME', ncid, dimid );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'SNCTOOLS:NC_ADD_RECS:inq_dimnameFailed', ncerr );
	end

	if dimid == -1
		error_id = 'SNCTOOLS:NC_ADD_RECS:noUnlimitedDimension';
		snc_error ( error_id, '%s is missing an unlimited dimension, %s requires it', ncfile, mfilename );
	end


else
	
	dimname = varargin{1};
	[dimid, status] = mexnc ( 'inq_dimid', ncid, dimname );
	if status ~= 0
		mexnc('close',ncid);
		ncerr = mexnc ( 'strerror', status );
		error_id = 'SNCTOOLS:NC_ADD_RECS:inq_dimidFailed';
		snc_error ( 'SNCTOOLS:NC_ADD_RECS:OPEN', ncerr );
	end
	
end
	
[dimlen, status] = mexnc ( 'INQ_DIMLEN', ncid, dimid );
if status ~= 0
	mexnc('close',ncid);
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_ADD_RECS:inq_dimlenFailed', ncerr );
end

status = mexnc('close',ncid);
if status ~= 0 
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_ADD_RECS:closeFailed', ncerr );
end

return



