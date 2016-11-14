function test_nc_info ( ncfile )
% TEST_NC_INFO:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests should fail
% test 1:  no input arguments, should fail
% test 2:  too many inputs
% test 3:  1 input, not a netcdf file
%
% 2nd set of tests should succeed
% test 4:  empty netcdf file
% test 5:  netcdf file has dimensions, but no variables.
% test 6:  netcdf file has singleton variables, but no dimensions.
% test 7:  netcdf file has global attributes, but no variables or dimensions
% test 8:  netcdf file with dimensions, variables, both unlimited variables 
%          and fixed variables, and global attributes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_info.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_INFO:  starting test suite...\n' );
test_001 ( ncfile );
test_002 ( ncfile );
test_003 ( ncfile );
test_004 ( ncfile );
test_005 ( ncfile );
test_006 ( ncfile );
test_007 ( ncfile );
test_008 ( ncfile );
return





function test_001 ( ncfile )
try
	nc = nc_info;
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_002 ( ncfile )
try
	nc = nc_info ( ncfile, 'blah' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_003 (ncfile)
fprintf ( 1, '    skipping test 003\n' );
%Ttry
%T	nc = nc_info ( '/dev/null' );
%T	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
%T	error ( msg );
%Tend
%fprintf ( 1, '%s:  passed.\n', testid );
return







function test_004 ( ncfile )

create_empty_file ( ncfile );
nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	msg = sprintf ( '%s:  :  Filename was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dimension ) ~= 0 )
	msg = sprintf ( '%s:  :  Dimension was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dataset ) ~= 0 )
	msg = sprintf ( '%s:  :  Dataset was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Attribute ) ~= 0 )
	msg = sprintf ( '%s:  :  Attribute was wrong.\n', mfilename  );
	error ( msg );
end
return









function test_005 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );
nc_add_dimension ( ncfile, 'x', 2 );
nc_add_dimension ( ncfile, 'y', 6 );

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	msg = sprintf ( '%s:  :  Filename was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dimension ) ~= 3 )
	msg = sprintf ( '%s:  :  Dimension was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dataset ) ~= 0 )
	msg = sprintf ( '%s:  :  Dataset was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Attribute ) ~= 0 )
	msg = sprintf ( '%s:  :  Attribute was wrong.\n', mfilename  );
	error ( msg );
end
return










function test_006 ( ncfile )


create_empty_file ( ncfile );

clear varstruct;
varstruct.Name = 'y';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	msg = sprintf ( '%s:  :  Filename was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dimension ) ~= 0 )
	msg = sprintf ( '%s:  :  Dimension was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dataset ) ~= 1 )
	msg = sprintf ( '%s:  :  Dataset was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Attribute ) ~= 0 )
	msg = sprintf ( '%s:  :  Attribute was wrong.\n', mfilename  );
	error ( msg );
end
return







function test_007 ( ncfile )

create_empty_file ( ncfile );

clear varstruct;
nc_attput ( ncfile, nc_global, 'test1', 'this' );
nc_attput ( ncfile, nc_global, 'test2', 'is' );
nc_attput ( ncfile, nc_global, 'test3', 'a' );
nc_attput ( ncfile, nc_global, 'test4', 'test' );

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	msg = sprintf ( '%s:  :  Filename was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dimension ) ~= 0 )
	msg = sprintf ( '%s:  :  Dimension was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dataset ) ~= 0 )
	msg = sprintf ( '%s:  :  Dataset was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Attribute ) ~= 4 )
	msg = sprintf ( '%s:  :  Attribute was wrong.\n', mfilename  );
	error ( msg );
end
return









function test_008 ( ncfile )

create_empty_file ( ncfile );

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
nc_attput ( ncfile, nc_global, 'test1', 'this' );
nc_attput ( ncfile, nc_global, 'test2', 'is' );
nc_attput ( ncfile, nc_global, 'test3', 'a' );
nc_attput ( ncfile, nc_global, 'test4', 'test' );

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	msg = sprintf ( '%s:  :  Filename was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dimension ) ~= 3 )
	msg = sprintf ( '%s:  :  Dimension was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Dataset ) ~= 7 )
	msg = sprintf ( '%s:  :  Dataset was wrong.\n', mfilename  );
	error ( msg );
end
if ( length ( nc.Attribute ) ~= 4 )
	msg = sprintf ( '%s:  :  Attribute was wrong.\n', mfilename  );
	error ( msg );
end
return







