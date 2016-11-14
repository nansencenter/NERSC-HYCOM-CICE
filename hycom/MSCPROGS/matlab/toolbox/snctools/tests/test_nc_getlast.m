function test_nc_getlast ( ncfile )
% TEST_NC_GETLAST:
%
% This first set of tests should all fail.
% Test 1:  No inputs.
% Test 2:  One input.
% Test 3:  Four inputs.
% Test 4:  1st input is not character.
% Test 5:  2nd input is not character.
% Test 6:  3rd input is not numeric.
% Test 7:  1st input is not a netcdf file.
% Test 8:  2nd input is not a netcdf variable.
% Test 9:  2nd input is a netcdf variable, but not unlimited.
% Test 10:  Non-positive "num_records"
% Test 11:  Time series variables have no data.
% Test 12:  Time series variables have data, but fewer than what was 
%           requested.
%
% This second set of tests should all succeed.
% Test 13:  Two inputs, should return the last record.
% Test 14:  Three valid inputs.
% Test 15:  Three valid inputs, but get everything.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_getlast.m 1806 2006-09-13 18:28:50Z johnevans007 $
% $LastChangedDate: 2006-09-13 14:28:50 -0400 (Wed, 13 Sep 2006) $
% $LastChangedRevision: 1806 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf ( 1, 'NC_GETLAST:  starting test suite...\n' );
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

return




function test_001 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast;
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return







function test_002 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nc_getlast ( ncfile );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_003 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( ncfile, 't1', 3, 4 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_004 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( 0, 't1' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_005 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( ncfile, 0 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_006 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( ncfile, 't1', 'a' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_007 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

fprintf ( 1, '    skipping test 007\n' );
%try
%	nb = nc_getlast ( '/dev/null', 't1', 1 );
%	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
%	error ( msg );
%end
%fprintf ( 1, '%s:  passed.\n', testid );
return




function test_008 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );


try
	nb = nc_getlast ( ncfile, 't4', 1 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_009 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( ncfile, 'x', 1 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_010 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

try
	nb = nc_getlast ( ncfile, 't1', 0 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return






function test_011 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );


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
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );




try
	nb = nc_getlast ( ncfile, 't1', 1 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );

return




function test_012 ( ncfile )


create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

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
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );


%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );


try
	nb = nc_getlast ( ncfile, 't1', 12 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return







function test_013 ( ncfile )
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

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
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );


%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );

v = nc_getlast ( ncfile, 't1' );
if ( length(v) ~= 1 )
	msg = sprintf ( '%s:  : return value length was wrong.\n', mfilename  );
	error ( msg );
end
return




function test_014 ( ncfile )
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

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
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );


%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );

v = nc_getlast ( ncfile, 't1', 7 );
if ( length(v) ~= 7 )
	msg = sprintf ( '%s:  : return value length was wrong.\n', mfilename  );
	error ( msg );
end
return



function test_015 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );

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
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'ocean_time' };
nc_addvar ( ncfile, varstruct );


%
% write ten records
x = [0:9]';
b.ocean_time = x;
b.t1 = x;
b.t2 = 1./(1+x);
b.t3 = x.^2;
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );

v = nc_getlast ( ncfile, 't1', 10 );
if ( length(v) ~= 10 )
	msg = sprintf ( '%s:  : return value length was wrong.\n', mfilename  );
	error ( msg );
end
return


