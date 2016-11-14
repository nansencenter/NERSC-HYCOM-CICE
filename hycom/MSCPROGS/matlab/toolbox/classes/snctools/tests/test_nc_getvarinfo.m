function test_nc_getvarinfo ( ncfile )
% TEST_NC_GETVARINFO:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests should fail
% test 1:  no input arguments, should fail
% test 2:  one input
% test 3:  too many inputs
% test 4:  2 inputs, 1st is not a netcdf file
% test 5:  2 inputs, 2nd is not a netcdf variable
% test 6:  2 inputs, 1st is character, 2nd is numeric
% test 7:  2 inputs, 2nd is character, 1st is numeric
%
% 2nd set of tests should succeed
% test 8:  empty netcdf variable with no data or attributes
% test 9:  empty netcdf variable (unlimited) with no data or attributes
% test 11:  1d netcdf variable with data and attributes
% test 12:  2d netcdf variable with data and attributes
% test 13:  character variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_getvarinfo.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 'NC_GETVARINFO:  starting test suite...\n' );

create_empty_file ( ncfile );

test_01 ( ncfile );
test_02 ( ncfile );
test_03 ( ncfile );
test_04 ( ncfile );
test_05 ( ncfile );
test_06 ( ncfile );
test_07 ( ncfile );

populate_ncfile ( ncfile );

test_08 ( ncfile );



%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );

test_09 ( ncfile );
test_11 ( ncfile );
test_12 ( ncfile );
test_13 ( ncfile );



return











function test_01 ( ncfile )

try
	nb = nc_getvarinfo;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end

return






function test_02 ( ncfile )

try
	nb = nc_getvarinfo ( ncfile );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end

return





function test_03 ( ncfile )
try
	nb = nc_getvarinfo ( ncfile, 't1' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
catch
	;
end

return





function test_04 ( ncfile )

try
	nb = nc_getvarinfo ( 'iamnotarealfilenoreally', 't1' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end

return






function populate_ncfile ( ncfile );

%
% make all the variable definitions.
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );
nc_add_dimension ( ncfile, 'y', 6 );

clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 'ocean_time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't1';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't2';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
varstruct.Attribute(1).Name = 'test_att';
varstruct.Attribute(1).Value = 'dud';
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 'y';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 'z';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'y', 'x' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 'c';
varstruct.Nctype = 'char';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );

return











function test_05 ( ncfile )

try
	nb = nc_getvarinfo ( ncfile, 't5' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return









function test_06 ( ncfile )
try
	nb = nc_getvarinfo ( ncfile, 0 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return




function test_07 ( ncfile )
try
	nb = nc_getvarinfo ( 0, 't1' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return




function test_08 ( ncfile )
v = nc_getvarinfo ( ncfile, 'x' );

if ~strcmp(v.Name, 'x' )
	msg = sprintf ( '%s:  Name was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Nctype~=6 )
	msg = sprintf ( '%s:  Nctype was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Unlimited~=0 )
	msg = sprintf ( '%s:  Unlimited was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Dimension)~=1 )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if ( ~strcmp(v.Dimension{1},'x') )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Size~=2 )
	msg = sprintf ( '%s:  Size was not correct.\n', mfilename  );
	error ( msg );
end
if (numel(v.Size)~=1 )
	msg = sprintf ( '%s:  Rank was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Attribute)~=0 )
	msg = sprintf ( '%s:  Attribute was not correct.\n', mfilename  );
	error ( msg );
end

return





function test_09 ( ncfile )

v = nc_getvarinfo ( ncfile, 't1' );

if ~strcmp(v.Name, 't1' )
	msg = sprintf ( '%s:  Name was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Nctype~=6 )
	msg = sprintf ( '%s:  Nctype was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Unlimited~=1 )
	msg = sprintf ( '%s:  Unlimited was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Dimension)~=1 )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Size~=10 )
	msg = sprintf ( '%s:  Size was not correct.\n', mfilename  );
	error ( msg );
end
if (numel(v.Size)~=1 )
	msg = sprintf ( '%s:  Rank was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Attribute)~=0 )
	msg = sprintf ( '%s:  Attribute was not correct.\n', mfilename  );
	error ( msg );
end

return







function test_11 ( ncfile )

v = nc_getvarinfo ( ncfile, 't2' );

if ~strcmp(v.Name, 't2' )
	msg = sprintf ( '%s:  Name was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Nctype~=6 )
	msg = sprintf ( '%s:  Nctype was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Unlimited~=1 )
	msg = sprintf ( '%s:  Unlimited was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Dimension)~=1 )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Size~=10 )
	msg = sprintf ( '%s:  Size was not correct.\n', mfilename  );
	error ( msg );
end
if (numel(v.Size)~=1 )
	msg = sprintf ( '%s:  Rank was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Attribute)~=1 )
	msg = sprintf ( '%s:  Attribute was not correct.\n', mfilename  );
	error ( msg );
end

return



function test_12 ( ncfile )

v = nc_getvarinfo ( ncfile, 'z' );

if ~strcmp(v.Name, 'z' )
	msg = sprintf ( '%s:  Name was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Nctype~=6 )
	msg = sprintf ( '%s:  Nctype was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Unlimited~=0 )
	msg = sprintf ( '%s:  Unlimited was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Dimension)~=2 )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if ( any(v.Size - [6 2]) )
	msg = sprintf ( '%s:  Size was not correct.\n', mfilename  );
	error ( msg );
end
if (numel(v.Size)~=2 )
	msg = sprintf ( '%s:  Rank was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Attribute)~=0 )
	msg = sprintf ( '%s:  Attribute was not correct.\n', mfilename  );
	error ( msg );
end

return














function test_13 ( ncfile )

v = nc_getvarinfo ( ncfile, 'c' );

if ~strcmp(v.Name, 'c' )
	msg = sprintf ( '%s:  Name was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Nctype~=2 )
	msg = sprintf ( '%s:  Nctype was not correct.\n', mfilename  );
	error ( msg );
end
if (v.Unlimited~=1 )
	msg = sprintf ( '%s:  Unlimited was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Dimension)~=1 )
	msg = sprintf ( '%s:  Dimension was not correct.\n', mfilename  );
	error ( msg );
end
if ( any(v.Size - [10]) )
	msg = sprintf ( '%s:  Size was not correct.\n', mfilename  );
	error ( msg );
end
if (numel(v.Size)~=1 )
	msg = sprintf ( '%s:  Rank was not correct.\n', mfilename  );
	error ( msg );
end
if (length(v.Attribute)~=0 )
	msg = sprintf ( '%s:  Attribute was not correct.\n', mfilename  );
	error ( msg );
end

return







