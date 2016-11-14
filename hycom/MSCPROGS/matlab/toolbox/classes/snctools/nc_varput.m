function nc_varput( ncfile, varname, data, varargin )
% NC_VARPUT:  Writes data into a netCDF file.
%
% NC_VARPUT(NCFILE,VARNAME,DATA) writes the matlab variable DATA to
% the variable VARNAME in the netCDF file NCFILE.  The main requirement
% here is that DATA have the same dimensions as the netCDF variable.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT) writes DATA contiguously, 
% starting at the zero-based index START and with extents given by
% COUNT.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT,STRIDE) writes DATA  
% starting at the zero-based index START with extents given by
% COUNT, but this time with strides given by STRIDE.  If STRIDE is not
% given, then it is assumes that all data is contiguous.
%
% EXAMPLES:
%    Suppose you have a netcdf variable called 'x' of size 6x4.  If you 
%    have an array of data called 'mydata' that is 6x4, then you can 
%    write to the entire variable with 
% 
%        >> nc_varput ( 'foo.nc', 'x', mydata );
%
%    If you wish to only write to the first 2 rows and three columns,
%    you could do the following
%
%        >> subdata = mydata(1:2,1:3);
%        >> nc_varput ( 'foo.nc', 'x', subdata, [0 0], [2 3] );
%
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_varput.m 2364 2007-10-08 20:28:39Z johnevans007 $
% $LastChangedDate: 2007-10-08 16:28:39 -0400 (Mon, 08 Oct 2007) $
% $LastChangedRevision: 2364 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


snc_nargchk(3,6,nargin);
snc_nargoutchk(0,0,nargout);


[start, count, stride] = parse_and_validate_args(ncfile,varname,varargin{:});



[ncid, status] = mexnc('open', ncfile, nc_write_mode);
if (status ~= 0)
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:OPEN', ncerr );
end




%
% check to see if the variable already exists.  
[varid, status] = mexnc('INQ_VARID', ncid, varname );
if ( status ~= 0 )
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VARID', ncerr );
end


[dud,var_type,nvdims,var_dim,dud, status]=mexnc('INQ_VAR',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VAR', ncerr );
end



v = nc_getvarinfo ( ncfile, varname );
nc_count = v.Size;


[start, count] = validate_indexing (ncid,nvdims,data,start,count,stride);

%
% check that the length of the start argument matches the rank of the variable.
if length(start) ~= length(nc_count)
    mexnc ( 'close', ncid );
    fmt = 'Length of START index (%d) does not make sense with a variable rank of %d.\n';
    msg = sprintf ( fmt, length(start), length(nc_count) );
    snc_error ( 'SNCTOOLS:NC_VARPUT:badIndexing', msg );
end



%
% Figure out which write routine we will use.  If the target variable is a singleton, then we must use
% VARPUT1.  If a stride was given, we must use VARPUTG.  Otherwise just use VARPUT.
if nvdims == 0
    write_op = 'put_var1';
elseif nargin == 3
    write_op = 'put_var';
elseif nargin == 5
    write_op = 'put_vara';
elseif nargin == 6
    write_op = 'put_vars';
else
    error ( 'unhandled write op.  How did we come to this??\n' );
end




data = handle_fill_value ( ncid, varid, data );
data = handle_scaling(ncid,varid,data);

write_the_data(ncid,varid,start,count,stride,write_op,data);


status = mexnc ( 'close', ncid );
if ( status ~= 0 )
    snc_error ( 'SNCTOOLS:nc_varput:close', mexnc('STRERROR',status));
end


return




function [start, count, stride] = parse_and_validate_args(ncfile,varname,varargin)

%
% Set up default outputs.
start = [];
count = [];
stride = [];


switch length(varargin)
case 2
    start = varargin{1};
    count = varargin{2};

case 3
    start = varargin{1};
    count = varargin{2};
    stride = varargin{3};

end



%
% Error checking on the inputs.
if ~ischar(ncfile)
    snc_error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the filename must be character.' );
end
if ~ischar(varname)
    snc_error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the variable name must be character.' );
end

if ~isnumeric ( start )
    snc_error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''start'' argument must be numeric.' );
end
if ~isnumeric ( count )
    snc_error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''count'' argument must be numeric.' );
end
if ~isnumeric ( stride )
    snc_error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''stride'' argument must be numeric.' );
end


return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [start, count] = validate_indexing(ncid,nvdims,data,start,count,stride)
% Check that any given start, count, and stride arguments actually make sense
% for this variable.  


%
% Singletons are a special case.  We need to set the start and count carefully.
if nvdims == 0

    if isempty(start) && isempty(count) && isempty(stride)

        %
        % This is the case of "nc_varput ( file, var, single_datum );"
        start = 0;
        count = 1;

    elseif ~isempty(start) && ~isempty(count) && ~isempty(stride)
        mexnc ( 'close', ncid );
        err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
        err_msg = 'Strides make no sense for a singleton variable.';
        snc_error ( err_id, err_msg );
    end

	return;

end

% If START and COUNT not given, and if not a singleton variable, then START is [0,..] and COUNT is 
% the size of the data.  
if isempty(start) && isempty(count) && ( nvdims > 0 )
    start = zeros(1,nvdims);
    count = zeros(1,nvdims);
	for j = 1:nvdims
		count(j) = size(data,j);
	end
end


%
% Check that the start, count, and stride arguments have the same length.
if ( numel(start) ~= numel(count) )
    mexnc ( 'close', ncid );
    err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
    err_msg = 'START and COUNT arguments must have the same length.';
    snc_error ( err_id, err_msg );
end
if ( ~isempty(stride) && (length(start) ~= length(stride)) )
    mexnc ( 'close', ncid );
    err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
    err_msg = 'START, COUNT, and STRIDE arguments must have the same length.';
    snc_error ( err_id, err_msg );
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_scaling(ncid,varid,data)
% HANDLE_MEX_SCALING
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.
%

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' );
if ( status == 0 )
    have_scale_factor = 1;
else
    have_scale_factor = 0;
end
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'add_offset' );
if ( status == 0 )
    have_add_offset = 1;
else
    have_add_offset = 0;
end

%
% Return early if we don't have either one.
if ~(have_scale_factor || have_add_offset)
    return;
end

scale_factor = 1.0;
add_offset = 0.0;


if have_scale_factor
	[scale_factor, status] = mexnc ( 'get_att_double', ncid, varid, 'scale_factor' );
	if ( status ~= 0 )
	    mexnc ( 'close', ncid );
	    ncerr = mexnc('strerror', status);
	    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', ncerr );
	end
end

if have_add_offset
	[add_offset, status] = mexnc ( 'get_att_double', ncid, varid, 'add_offset' );
	if ( status ~= 0 )
	    mexnc ( 'close', ncid );
	    ncerr = mexnc('strerror', status);
	    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', ncerr );
	end
end

[var_type,status]=mexnc('INQ_VARTYPE',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VARTYPE', ncerr );
end

data = (double(data) - add_offset) / scale_factor;

%
% When scaling to an integer, we should add 0.5 to the data.  Otherwise
% there is a tiny loss in precision, e.g. 82.7 should round to 83, not 
% .
switch var_type
	case { nc_int, nc_short, nc_byte, nc_char }
	    data = round(data);
end


return






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_fill_value(ncid,varid,data)

%
% Handle the fill value.  We do this by changing any NaNs into
% the _FillValue.  That way the netcdf library will recognize it.
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' );
if ( status == 0 )

    switch ( class(data) )
    case 'double'
		funcstr = 'get_att_double';
    case 'single'
		funcstr = 'get_att_float';
    case 'int32'
		funcstr = 'get_att_int';
    case 'int16'
		funcstr = 'get_att_short';
    case 'int8'
		funcstr = 'get_att_schar';
    case 'uint8'
		funcstr = 'get_att_uchar';
    case 'char'
		funcstr = 'get_att_text';
    otherwise
	    mexnc ( 'close', ncid );
	    msg = sprintf ( 'Unhandled datatype for fill value, ''%s''.', class(data) );
	    error ( 'SNCTOOLS:NC_VARPUT:unhandledDatatype', msg );
    end

    [fill_value, status] = mexnc(funcstr,ncid,varid,'_FillValue' );
    if ( status ~= 0 )
	    mexnc ( 'close', ncid );
	    ncerr = mexnc('strerror', status);
		err_id = [ 'SNCTOOLS:NC_VARPUT:MEXNC:' funcstr ];
	    snc_error ( err_id, ncerr );
    end


    data(isnan(data)) = fill_value;

end

    













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_the_data(ncid,varid,start,count,stride,write_op,data)

%
% write the data
pdata = permute(data, fliplr( 1:ndims(data) ));
switch ( write_op )

    case 'put_var1'
        switch ( class(data) )
        case 'double'
			funcstr = 'put_var1_double';
        case 'single'
			funcstr = 'put_var1_float';
        case 'int32'
			funcstr = 'put_var1_int';
        case 'int16'
			funcstr = 'put_var1_short';
        case 'int8'
			funcstr = 'put_var1_schar';
        case 'uint8'
			funcstr = 'put_var1_uchar';
        case 'char'
			funcstr = 'put_var1_text';
        otherwise
            mexnc('close',ncid);
            msg = sprintf ( 'unhandled data class %s\n', class(pdata) );
            snc_error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', msg );
        end
        status = mexnc (funcstr, ncid, varid, start, pdata );

    case 'put_var'
        switch ( class(data) )
        case 'double'
			funcstr = 'put_var_double';
        case 'single'
			funcstr = 'put_var_float';
        case 'int32'
			funcstr = 'put_var_int';
        case 'int16'
			funcstr = 'put_var_short';
        case 'int8'
			funcstr = 'put_var_schar';
        case 'uint8'
			funcstr = 'put_var_uchar';
        case 'char'
			funcstr = 'put_var_text';
        otherwise
            mexnc('close',ncid);
            msg = sprintf ( 'unhandled data class %s\n', class(pdata) );
            snc_error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', msg );
        end
        status = mexnc (funcstr, ncid, varid, pdata );
    
    case 'put_vara'
        switch ( class(data) )
        case 'double'
			funcstr = 'put_vara_double';
        case 'single'
			funcstr = 'put_vara_float';
        case 'int32'
			funcstr = 'put_vara_int';
        case 'int16'
			funcstr = 'put_vara_short';
        case 'int8'
			funcstr = 'put_vara_schar';
        case 'uint8'
			funcstr = 'put_vara_uchar';
        case 'char'
			funcstr = 'put_vara_text';
        otherwise
            mexnc('close',ncid);
            msg = sprintf ( 'unhandled data class %s\n', class(pdata) );
            snc_error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', msg );
        end
        status = mexnc (funcstr, ncid, varid, start, count, pdata );

    case 'put_vars'
        switch ( class(data) )
        case 'double'
			funcstr = 'put_vars_double';
        case 'single'
			funcstr = 'put_vars_float';
        case 'int32'
			funcstr = 'put_vars_int';
        case 'int16'
			funcstr = 'put_vars_short';
        case 'int8'
			funcstr = 'put_vars_schar';
        case 'uint8'
			funcstr = 'put_vars_uchar';
        case 'char'
			funcstr = 'put_vars_text';
        otherwise
            mexnc('close',ncid);
            msg = sprintf ( 'unhandled data class %s\n', class(pdata) );
            snc_error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', msg );
        end
        status = mexnc (funcstr, ncid, varid, start, count, stride, pdata );

    otherwise 
        mexnc ( 'close', ncid );
        msg = sprintf ( 'unknown write operation''%s''.\n', write_op );
        snc_error ( 'SNCTOOLS:NC_VARPUT:unhandledWriteOp', msg );


end

if ( status ~= 0 )
    mexnc ( 'close', ncid );
    ncerr = mexnc ( 'strerror', status );
    msg = sprintf ( 'write operation ''%s'' failed with error ''%s''.\n', ...
					write_op, ncerr );
    error ( msg );
end

return
