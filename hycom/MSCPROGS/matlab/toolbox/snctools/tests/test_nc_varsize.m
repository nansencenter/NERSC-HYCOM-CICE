function test_nc_varsize ( ncfile )
% TEST_NC_VARSIZE:
%
% Depends upon nc_add_dimension, nc_addvar
%
% 1st set of tests, routine should fail
% test 1:  no input arguments
% test 2:  1 input
% test 3:  too many inputs
% test 4:  inputs are not all character
% test 5:  not a netcdf file
% test 6:  empty netcdf file
% test 7:  given variable is not present
%
% 2nd set of tests, routine should succeed
% test 8:  given singleton variable is present
% test 9:  given 1D variable is present
% test 10:  given 1D-unlimited-but-empty variable is present
% test 11:  given 2D variable is present

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_varsize.m 1806 2006-09-13 18:28:50Z johnevans007 $
% $LastChangedDate: 2006-09-13 14:28:50 -0400 (Wed, 13 Sep 2006) $
% $LastChangedRevision: 1806 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf ( 1, 'NC_VARSIZE:  starting test suite...\n' );

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

return



function test_001 ( ncfile )

try
	v = nc_varsize;
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return











function test_002 ( ncfile )

try
	v = nc_varsize ( ncfile );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return










function test_003 ( ncfile )

try
	v = nc_varsize ( ncfile, 'x', 'y' );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return










function test_004 ( ncfile )

try
	v = nc_varsize ( ncfile, 1 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return











function test_005 ( ncfile )

% test 5:  not a netcdf file
fprintf ( 1, '    skipping test 005...\n' );
%try
%	v = nc_varsize ( '/dev/null', 't' );
%	msg = sprintf ( '%s:  %s succeeded when it should have failed.\n', mfilename, testid );
%	error ( msg );
%end
%fprintf ( 1, '%s:  passed.\n', testid );
return














function test_006 ( ncfile )

create_empty_file ( ncfile );
try
	v = nc_varsize ( ncfile, 't' );
	msg = sprintf ( '%s:  succceeded when it should have failed.\n', mfilename );
	error ( msg );
end
return










function test_007 ( ncfile )

create_empty_file ( ncfile );
try
	v = nc_varsize ( ncfile, 'x' );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return











function test_008 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 's';
varstruct.Nctype = 'double';
varstruct.Dimension = [];
nc_addvar ( ncfile, varstruct );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

varsize = nc_varsize ( ncfile, 's' );
if ( varsize ~= 1 )
	error ( '%s:  varsize was not right.\n', mfilename );
end
return









function test_009 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 's';
varstruct.Nctype = 'double';
varstruct.Dimension = { 's' };
nc_addvar ( ncfile, varstruct );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

varsize = nc_varsize ( ncfile, 's' );
if ( varsize ~= 5 )
	error ( '%s:  varsize was not right.\n', mfilename );
end
return











function test_010 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 't';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't', 's' };
nc_addvar ( ncfile, varstruct );
clear varstruct;
varstruct.Name = 'w';
varstruct.Nctype = 'double';
varstruct.Dimension = { 't' };
nc_addvar ( ncfile, varstruct );

varsize = nc_varsize ( ncfile, 't' );
if ( varsize(1) ~= 0 ) & ( varsize(2) ~= 5 )
	error ( '%s:  varsize was not right.\n', mfilename );
end
return










function test_011 ( ncfile )


create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 's', 5 );
nc_add_dimension ( ncfile, 't', 7 );
clear varstruct;
varstruct.Name = 'st';
varstruct.Nctype = 'double';
varstruct.Dimension = { 's', 't' };
nc_addvar ( ncfile, varstruct );

varsize = nc_varsize ( ncfile, 'st' );
if ( varsize(1) ~= 5 ) & ( varsize(2) ~= 7 )
	error ( '%s:  varsize was not right.\n', mfilename );
end
return









