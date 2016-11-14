function test_nc_varrename ( ncfile )
% TEST_NC_ISVAR:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests, routine should fail
% test 1:  no input arguments
% test 2:  1 input
% test 3:  2 input
% test 4:  too many inputs
% test 5:  inputs are not all character
% test 6:  not a netcdf file
% test 7:  empty netcdf file
% test 8:  given variable is not present
%
% 2nd set of tests, routine should succeed
% test 9:  given variable is present

%
% 3rd set should fail
% test 10:  given variable is present, but another exists with the same name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_varrename.m 1806 2006-09-13 18:28:50Z johnevans007 $
% $LastChangedDate: 2006-09-13 14:28:50 -0400 (Wed, 13 Sep 2006) $
% $LastChangedRevision: 1806 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_VARRENAME:  starting test suite.\n' );

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
	nc_varrename;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return












function test_002 ( ncfile )

try
	nc_varrename ( ncfile );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return








function test_003 ( ncfile )

try
	nc_varrename ( ncfile, 'x' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return









function test_004 ( ncfile )

try
	nc_varrename ( ncfile, 'blah', 'blah2', 'blah3' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return







function test_005 ( ncfile )

%
% Ok, now we'll create the test file
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );





% test 5:  inputs are not all character
try
	nc_varrename ( ncfile, 'x', 1 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end

return










function test_006 ( ncfile )
% test 6:  not a netcdf file
fprintf ( 1, '    skipping test 006.\n' );
%try
%	nc_varrename ( '/dev/null', 't', 'u' );
%	msg = sprintf ( '%s:  %s succeeded when it should have failed.\n', mfilename, testid );
%	error ( msg );
%end
%fprintf ( 1, '%s:  passed.\n', testid );
return












function test_007 ( ncfile )


create_empty_file ( ncfile );
try
	nc_varrename ( ncfile, 'x', 'y' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return








function test_008 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

try
	nc_varrename ( ncfile, 't2', 't3' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
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

nc_varrename ( ncfile, 't', 't2' );



v = nc_getvarinfo ( ncfile, 't2' );
if ~strcmp ( v.Name, 't2' )
	msg = sprintf ( '%s: rename did not seem to work.\n', mfilename );
	error ( msg );
end

return










function test_010 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );
varstruct.Name = 't2';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

try
	nc_varrename ( ncfile, 't', 't2' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end



v = nc_getvarinfo ( ncfile, 't2' );
if ~strcmp ( v.Name, 't2' )
	msg = sprintf ( '%s:  rename did not seem to work.\n', mfilename );
	error ( msg );
end


return


