function test_nc_iscoordvar ( ncfile )
% TEST_NC_ISCOORDVAR:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests should fail
% test 1:  no input arguments
% test 2:  1 input
% test 3:  too many inputs
% test 4:  both inputs are not character
% test 5:  not a netcdf file
% test 6:  empty netcdf file
% test 7:  netcdf file has dimensions, but no variables.
% test 8:  given variable is not present  
% test 9:  given variable's dimension is not of the same name
% test 10:  given variable has more than one dimension
%
% 2nd set of tests should succeed
% test 11:  netcdf file has singleton variable, but no dimensions.
% test 12:  given variable has one dimension of the same name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_iscoordvar.m 1806 2006-09-13 18:28:50Z johnevans007 $
% $LastChangedDate: 2006-09-13 14:28:50 -0400 (Wed, 13 Sep 2006) $
% $LastChangedRevision: 1806 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf ( 1, 'NC_ISCOORDVAR:  starting test suite...\n' );
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

return









function test_001 ( ncfile )
try
	nc = nc_iscoordvar;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end




function test_002 ( ncfile )

try
	nc = nc_iscoordvar ( ncfile );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_003 ( ncfile )

try
	nc = nc_iscoordvar ( ncfile, 'blah', 'blah2' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
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





% test 4:  both inputs are not character
testid = 'Test 4';
fprintf ( 1, '    skipping test 004.\n' );
%try
%	nc = nc_iscoordvar ( '/dev/null', 5 );
%	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
%	error ( msg );
%end
%fprintf ( 1, '%s:  passed.\n', testid );
return






function test_005 ( ncfile )
% test 5:  not a netcdf file
fprintf ( 1, '    skipping test 005.\n' );
%try
%	nc = nc_iscoordvar ( '/dev/null', 't' );
%	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
%%	error ( msg );
%end
%fprintf ( 1, '%s:  passed.\n', testid );
return







function test_006 ( ncfile )

create_empty_file ( ncfile );
try
	nc = nc_iscoordvar ( ncfile, 't' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return










function test_007 (ncfile )
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
try
	nc = nc_iscoordvar ( ncfile, 't' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return







function test_008 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = [];

nc_addvar ( ncfile, varstruct );
try
	nc = nc_iscoordvar ( ncfile, 'y' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return








function test_009 ( ncfile )

% 2nd set of tests should succeed
% test 9:  given variable's dimension is not of the same name

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 's' };

nc_addvar ( ncfile, varstruct );
b = nc_iscoordvar ( ncfile, 't' );
if ( b ~= 0 )
	msg = sprintf ( '%s:    incorrect result.\n', mfilename  );
	error ( msg );
end
return






function test_010 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 7 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 's', 't' };
nc_addvar ( ncfile, varstruct );

b = nc_iscoordvar ( ncfile, 't' );
if ( b ~= 0 )
	msg = sprintf ( '%s:    incorrect result.\n', mfilename  );
	error ( msg );
end
return







function test_011 ( ncfile )

%
% test 11:  netcdf file has singleton variable, but no dimensions.
testid = 'Test 11';
create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );

b = nc_iscoordvar ( ncfile, 't' );
return






function test_012 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

b = nc_iscoordvar ( ncfile, 't' );

return









