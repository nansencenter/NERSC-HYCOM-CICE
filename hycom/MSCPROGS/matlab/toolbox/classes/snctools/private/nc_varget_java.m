function [values, the_var_size] = nc_varget_java ( ncfile, varname, start, count, stride )
% NC_VARGET_JAVA:  Java backend for nc_varget.
%
% See the help section for nc_varget.


import ucar.nc2.dods.*     % import opendap reader classes
import ucar.nc2.*          % have to import this (NetcdfFile) 
                           % as well for local reads.
                           % Now that's just brilliant.  






if snc_is_url ( ncfile )
	try  
		jncid = DODSNetcdfFile(ncfile);
	catch
		try
			jncid = NetcdfFile.open(ncfile);
		catch
			error ( 'SNCTOOLS:nc_varget_java:urlOpen', 'Could not open netCDF URL.' );
		end
	end
else
	try
		jncid = NetcdfFile.open(ncfile);
	catch
		error ( 'SNCTOOLS:nc_varget_java:localFileOpen', 'Could not open netCDF file.' );
	end
end




%
% Get the variable object
jvarid = jncid.findVariable(varname);
if isempty ( jvarid )
	msg = sprintf ('findVariable failed on variable ''%s'', file ''%s''.', varname, ncfile );
    snc_error ( 'SNCTOOLS:NC_VARGET:JAVA:noSuchVariable', msg );
end



theDataType = jvarid.getDataType();
theDataTypeString = char ( theDataType.toString() ) ;
theDimensions = jvarid.getDimensions();
num_var_dims = theDimensions.size();


%
% Check that the start, count, stride parameters have appropriate lengths.
% Otherwise we risk confusing the mex-file.
validate_index_vectors ( start, count, stride, num_var_dims );



varinfo = nc_getvarinfo ( ncfile, varname );
the_var_size = varinfo.Size;



[read_method, start, count] = determine_read_method_java(start,count,stride,varinfo);



%
% If the user had set non-positive numbers in "count", then
% we replace them with what we need to get the
% rest of the variable.
negs = find(count<0);
count(negs) = the_var_size(negs) - start(negs);






%
% Finally!  Read the freakin' data.
try
    switch ( read_method )
        case 'GET_VAR1'
    
            % Read a scalar
            values = nc_var1_get_java ( jvarid, theDataTypeString );
    
        case 'GET_VAR'
        
            % Read everything.
            values = nc_var_get_java ( jvarid, theDataTypeString );
    
        case 'GET_VARA'
        
            % Read a contiguous subset
            values = nc_vara_get_java ( jvarid, theDataTypeString, start, count );
    
        case 'GET_VARS'
        
            % Read a contiguous subset
            values = nc_vars_get_java ( jvarid, theDataTypeString, start, count, stride );
    
        otherwise
        
            snc_error ( 'SNCTOOLS_NC_VARGET:JAVA:badReadCase', 'Unhandled case, don''t know which read method to use.' );  
        
    end
    
catch
    close ( jncid );
    rethrow ( lasterror );
end
    
    

%
% If it's a 1D vector, make it a column vector.
% Apparently, there's no need to otherwise permute the data.
if length(the_var_size) == 1
    values = values(:);
end                                                                                   



values = handle_fill_value_java ( jvarid, theDataType, values );
values = handle_missing_value_java ( jvarid, theDataType, values );
values = handle_scaling_java ( jvarid, values );


%
% remove any singleton dimensions.
values = squeeze ( values );



close ( jncid );


return









function values = nc_var1_get_java ( jvarid, theDataTypeString )
% NC_VAR1_GET_JAVA:  reads a scalar.

switch ( theDataTypeString )

    case 'char'
        values = jvarid.read();
        values = char ( values.toString() );

    case { 'double', 'float', 'int', 'short', 'byte' }
        values = jvarid.readScalarDouble();

    otherwise
        msg = sprintf ('unhandled datatype ''%s''', theDataTypeString );;
        error ( msg );
    
end
    
    
return






function values = nc_var_get_java ( jvarid, theDataTypeString )
% NC_VAR_GET_JAVA:  reads the entire variable

values = jvarid.read();
switch ( theDataTypeString )
    case 'char'
        values = copyToNDJavaArray(values);
    case { 'double', 'float', 'int', 'short', 'byte' }
        values = copyToNDJavaArray(values);
        values = double(values);
    otherwise
        msg = sprintf ( 'unhandled datatype ''%s''', theDataTypeString );
        error ( msg );
    
end
return







    
    
function values = nc_vara_get_java ( jvarid, theDataTypeString, start, count )
% NC_VARA_GET_JAVA:  reads a contiguous subset

values = jvarid.read(start, count);
switch ( theDataTypeString )
    case 'char'
        values = copyToNDJavaArray(values);
    case { 'double', 'float', 'int', 'short', 'byte' }
        values = copyToNDJavaArray(values);
        values = double ( values );
    otherwise
        msg = sprintf ( 'unhandled datatype ''%s''', theDataTypeString );
        error (  msg );
    
end
return
    







function values = nc_vars_get_java ( jvarid, theDataTypeString, start, count, stride )
% NC_VARS_GET_JAVA:  reads a strided subset

% Have to use the method with the section selector.
% "1:2,10,:,1:100:10"
extent = start + count.*stride-1;
section_selector = '';
for j = 1:length(start)
    section_selector = sprintf ( '%s,%d:%d:%d', section_selector, start(j), extent(j), stride(j) );
end

%
% Get rid of the first comma.
section_selector(1) = [];

values = jvarid.read(section_selector);
switch ( theDataTypeString )
    case 'char'
        values = copyToNDJavaArray(values);
    case { 'double', 'float', 'int', 'short', 'byte' }
        values = copyToNDJavaArray(values);
        values = double ( values );
    otherwise
        msg = sprintf ( 'unhandled datatype ''%s''', theDataTypeString );
        error ( msg );
    
end
    
    
return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DETERMINE_READ_METHOD_JAVA
%    Determine the read method that we will instruct java to use in order to
%    properly read in the netCDF data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [read_method, start, count] = determine_read_method_java ( start, count, stride, varinfo )

%
% If a singleton, then use GET_VAR1.  This is only because some 
% opendap-enabled mexnc clients have trouble using GET_VAR on 
% singletons.  It is annoying to have to do this, but it works just as
% well.
if length(varinfo.Dimension) == 0
	read_method = 'GET_VAR1';
	start = 0;
	count = 1;
	return
end


if isempty(start) && isempty(count) && isempty(stride)
    read_method = 'GET_VAR';
elseif ~isempty(start) && ~isempty(count) && isempty(stride)
    read_method = 'GET_VARA';
elseif ~isempty(start) && ~isempty(count) && ~isempty(stride)
    read_method = 'GET_VARS';
else
    snc_error ( 'SNCTOOLS:NC_VARGET:JAVA:undeterminedReadCase', ...
            'Could not determine intended read method.' );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HANDLE_FILL_VALUE_JAVA
%     If there is a fill value, then replace such values with NaN.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_fill_value_java ( jvarid, var_type, values )

%
% Handle the fill value, if any.  Change those values into NaN.
fillvalue_att = jvarid.findAttribute ( '_FillValue' );
if ~isempty(fillvalue_att)
    switch ( char ( var_type.toString() ) )
    case 'char'
        %
        % For now, do nothing.  Does a fill value even make sense with char data?
        % If it does, please tell me so.
        %

    case { 'double', 'float', 'long', 'short', 'byte' }
        fill_value = fillvalue_att.getNumericValue().doubleValue();
        values(values==fill_value) = NaN;

    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HANDLE_MISSING_VALUE_JAVA
%     If there is a missing value, then replace such values with NaN.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_missing_value_java ( jvarid, theDataType, values )

%
% If there is a fill value attribute, then that had precedence.  Do nothing.
fvatt = jvarid.findAttribute ( '_FillValue' );
if ~isempty(fvatt)
	return
end

%
% Handle the missing value, if any.  Change those values into NaN.
missing_value_att = jvarid.findAttribute ( 'missing_value' );
if ~isempty(missing_value_att)
    switch ( char ( theDataType.toString() ) )
    case 'char'

        %
        % For now, do nothing.  Does a fill value even make sense with 
        % char data?  Matlab doesn't allow for NaNs in character arrays.

    case { 'double', 'float', 'long', 'short', 'byte' }
        missing_value = missing_value_att.getNumericValue().doubleValue();
        values(values==missing_value) = NaN;

    end
end

return









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HANDLE_SCALING_JAVA
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = handle_scaling_java ( jvarid, values )

%
% Handle the scale factor and add_offsets. 
scale_factor_att = jvarid.findAttribute ( 'scale_factor' );
add_offset_att = jvarid.findAttribute ( 'add_offset' );

%
% Return early if we don't have either one.
if isempty(scale_factor_att) && isempty(add_offset_att)
    return
end


if ~isempty(scale_factor_att)
    scale_factor = scale_factor_att.getNumericValue().doubleValue();
else
    scale_factor = 1.0;
end

if ~isempty(add_offset_att)
    add_offset = add_offset_att.getNumericValue().doubleValue();
else
    add_offset = 0.0;
end


values = values * scale_factor + add_offset;

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DETERMINE_VARSIZE_JAVA
%    Need to figure out just how big the variable is.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function the_var_size = determine_varsize_java ( jvarid, num_var_dims )
%
% If a singleton, then use GET_VAR1.  This is only because some 
% opendap-enabled mexnc clients have trouble using GET_VAR on 
% singletons.  It is annoying to have to do this, but it works just as
% well.
if num_var_dims == 0

    %
    % Must be a scalar variable.  
    the_var_size = 1;

else
    the_var_size = zeros(1,num_var_dims);
    for j=1:num_var_dims,
            the_var_size(j) = jvarid.getDimension(j-1).getLength();
    end
end

return






