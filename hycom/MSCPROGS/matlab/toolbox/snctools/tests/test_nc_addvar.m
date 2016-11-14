function test_nc_addvar ( ncfile )
% TEST_NC_ADDVAR
%
% Relies upon nc_getvarinfo.
%
% Need to have more tests
% Test 1:  no inputs
% test 2:  too many inputs
% test 3:  first input not a netcdf file
% test 4:  2nd input not character
% test 5:  3rd input not numeric
% test 6:  Add a double variable, no dimensions
% test 7:  Add a float variable
% test 8:  Add a int32 variable
% test 9:  Add a int16 variable
% test 10:  Add a byte variable
% test 11:  Add a char variable
% test 12:  Add a variable with a named fixed dimension.
% test 13:  Add a variable with a named unlimited dimension.
% test 14:  Add a variable with a named unlimited dimension and named fixed dimension.
% test 15:  Add a variable with attribute structures.
% test 17:  Add a variable with a numeric Nctype
% test 18:  Try to add a variable that is already there.
% test 19:  Have an illegal field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_addvar.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ADDVAR:  starting test suite...\n' );

test_001 ( ncfile );
test_002 ( ncfile );
test_003 ( ncfile );
test_004 ( ncfile );
test_005 ( ncfile );
test_006 ( ncfile );
test_007 ( ncfile );
test_008 ( ncfile );
test_009 ( ncfile );
test_010 ( ncfile );
test_011 ( ncfile );
test_012 ( ncfile );
test_013 ( ncfile );
test_014 ( ncfile );
test_015 ( ncfile );
%test_016 ( ncfile );
test_017 ( ncfile );
test_018 ( ncfile );
test_019 ( ncfile );

return













function test_001 ( ncfile );

try
	nc_addvar;
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return












function test_002 ( ncfile )
% test 2:  too many inputs
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'test_var';
try
	nc_addvar ( ncfile, varstruct );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
catch
	;
end
return














function test_003 ( ncfile )
% test 3:  first input not a netcdf file
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'test_var';
fprintf ( 1, '%s:  skipping test_003.\n' );
%try
%	nc_addvar ( '/dev/null', varstruct );
%	msg = sprintf ( '%s:  %s succeeded when it should have failed.\n', mfilename, testid );
%	error ( msg );
%end
return






function test_004 ( ncfile )
% test 4:  2nd input not character
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'test_var';
try
	nc_addvar ( ncfile, 5 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return













function test_005 ( ncfile )
% test 5:  No name provided in varstruct
create_empty_file ( ncfile );
clear varstruct;
varstruct = struct([]);
try
	nc_addvar ( ncfile, varstruct );
	error ( '%s:  %s succeeded when it should have failed.\n', mfilename );
end
return










function test_006 ( ncfile )
% test 6:  Add a double variable, no dimensions
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end
if ( v.Size ~= 1 ) 
	error ( '%s:  data size was wrong.\n', mfilename );
end
if ( ~isempty(v.Dimension) ) 
	error ( '%s:  dimensions were wrong.\n', mfilename );
end

return










function test_007 ( ncfile )

% test 7:  Add a float variable
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'float';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_FLOAT' )
	error ( '%s:  data type was wrong.\n', mfilename );
end

return






function test_008 ( ncfile )

% test 8:  Add a int32 variable
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'int';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_INT' )
	error ( '%s:  data type was wrong.\n', mfilename );
end

return






function test_009 ( ncfile )
% test 9:  Add a int16 variable
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'short';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_SHORT' )
	error ( '%s:  data type was wrong.\n', mfilename );
end

return









function test_010 ( ncfile )

% test 10:  Add a byte variable
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'byte';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_BYTE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end

return









function test_011 ( ncfile )

% test 11:  Add a char variable
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'char';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_CHAR' )
	error ( '%s:  data type was wrong.\n', mfilename );
end

return







function test_012 ( ncfile )

% test 12:  Add a variable with a named fixed dimension.
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 5 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end
if any(v.Size - 5)
	error ( '%s:  variable size was wrong.\n', mfilename );
end
if ( length(v.Dimension) ~= 1 ) 
	error ( '%s:  dimensions were wrong.\n', mfilename );
end
if ( ~strcmp(v.Dimension{1}, 'x' ) ) 
	error ( '%s:  dimensions were wrong.\n', mfilename );
end

return










function test_013 ( ncfile )

% test 13:  Add a variable with a named unlimited dimension.
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 0 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end
if (v.Size ~= 0)
	error ( '%s:  variable size was wrong.\n', mfilename );
end
if ( ~v.Unlimited)
	error ( '%s:  unlimited classifaction was wrong.\n', mfilename );
end

return












function test_014 ( ncfile )
% test 14:  Add a variable with a named unlimited dimension and named fixed dimension.
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 0 );
nc_add_dimension ( ncfile, 'y', 5 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x', 'y' };
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  %s:  data type was wrong.\n', mfilename, testid );
end
if (v.Size(1) ~= 0) & (v.Size(2) ~= 5 )
	error ( '%s:  %s:  variable size was wrong.\n', mfilename, testid );
end
if ( ~v.Unlimited )
	error ( '%s:  %s:  unlimited classifaction was wrong.\n', mfilename, testid );
end

return









function test_015 ( ncfile )
% test 15:  Add a variable with attribute structures.
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 0 );
nc_add_dimension ( ncfile, 'y', 5 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };
varstruct.Attribute(1).Name = 'test';
varstruct.Attribute(1).Value = 'blah';
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end
if (v.Size(1) ~= 0) & (v.Size(2) ~= 5 )
	error ( '%s:  variable size was wrong.\n', mfilename );
end
if ( ~v.Unlimited)
	error ( '%s:  unlimited classifaction was wrong.\n', mfilename );
end
if ( length(v.Attribute) ~= 1)
	error ( '%s:  number of attributes was wrong.\n', mfilename );
end

return









function test_017 ( ncfile )

% test 17:  Add a variable with a numeric Nctype
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 5 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 6;
varstruct.Dimension = { 'x' };
nc_addvar ( ncfile, varstruct );

v = nc_getvarinfo ( ncfile, 'x' );
if ~strcmp(nc_datatype_string(v.Nctype),'NC_DOUBLE' )
	error ( '%s:  data type was wrong.\n', mfilename );
end
if any(v.Size - 5)
	error ( '%s:  variable size was wrong.\n', mfilename );
end
if ( length(v.Dimension) ~= 1 ) 
	error ( '%s:  dimensions were wrong.\n', mfilename );
end
if ( ~strcmp(v.Dimension{1}, 'x' ) ) 
	error ( '%s:  dimensions were wrong.\n', mfilename );
end




return












function test_018 ( ncfile )


create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 5 );

clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 6;
varstruct.Dimension = { 'x' };
nc_addvar ( ncfile, varstruct );
try
	nc_addvar ( ncfile, varstruct );
	msg = sprintf ( '%s:  succeeded when it should have failed', mfilename );
end



return



function test_019 ( ncfile )


create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'x', 5 );

clear varstruct;
varstruct.Name = 'x';
varstruct.nnccttyyppee = { 'x' };
varstruct.Dimension = { 'x' };
try
	nc_addvar ( ncfile, varstruct );
	msg = sprintf ( '%s:  succeeded when it should have failed', mfilename );
end



return
















