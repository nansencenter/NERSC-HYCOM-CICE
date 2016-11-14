function test_nc_isvar ( ncfile )
% TEST_NC_ISVAR:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests, routine should fail
% test 1:  no input arguments
% test 2:  1 input
% test 3:  too many inputs
% test 4:  both inputs are not character
% test 5:  not a netcdf file
% test 6:  empty netcdf file
% test 7:  netcdf file has dimensions, but no variables.
%
% 2nd set of tests, routine should succeed
% test 8:  given variable is not present
% test 9:  given 1D variable is present

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_isvar.m 1806 2006-09-13 18:28:50Z johnevans007 $
% $LastChangedDate: 2006-09-13 14:28:50 -0400 (Wed, 13 Sep 2006) $
% $LastChangedRevision: 1806 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ISVAR:  starting test suite...\n' );

test_001 ( ncfile );
test_002 ( ncfile );
test_003 ( ncfile );
test_004 ( ncfile );
test_005 ( ncfile );
test_006 ( ncfile );
test_007 ( ncfile );
test_008 ( ncfile );
test_009 ( ncfile );

return





function test_001 ( ncfile )

try
	nc = nc_isvar;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return









function test_002 ( ncfile )

try
	nc = nc_isvar ( ncfile );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return










function test_003 ( ncfile )

try
	nc = nc_isvar ( ncfile, 'blah', 'blah2' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return









function test_004 ( ncfile )

%
% Ok, now we'll create the test file
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );





try
	nc = nc_isvar ( ncfile, 5 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return













function test_005 ( ncfile )

% test 5:  not a netcdf file
fprintf ( 1, '    skipping test 005.\n' );
%try
%	nc = nc_isvar ( '/dev/null', 't' );
%	msg = sprintf ( '%s:  %s succeeded when it should have failed.\n', mfilename, testid );
%	error ( msg );
%end
return










function test_006 ( ncfile )

create_empty_file ( ncfile );
yn = nc_isvar ( ncfile, 't' );
if ( yn == 1 )
	msg = sprintf ( '%s:  incorrectly classified.\n', mfilename );
	error ( msg );
end
return











function test_007 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
yn = nc_isvar ( ncfile, 't' );
if ( yn == 1 )
	msg = sprintf ( '%s:  incorrectly classified.\n', mfilename );
	error ( msg );
end
return













function test_008 ( ncfile )


create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 0 );
nc_add_dimension ( ncfile, 'x', 10 );
clear varstruct;
varstruct.Name = 'x';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };

nc_addvar ( ncfile, varstruct );

b = nc_isvar ( ncfile, 'y' );
if ( b ~= 0 )
	msg = sprintf ( '%s:  incorrect result.\n', mfilename );
	error ( msg );
end
return











function test_009 ( ncfile )



create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

b = nc_isvar ( ncfile, 't' );
if ( b ~= 1 )
	msg = sprintf ( '%s:   incorrect result.\n', mfilename );
	error ( msg );
end
return
