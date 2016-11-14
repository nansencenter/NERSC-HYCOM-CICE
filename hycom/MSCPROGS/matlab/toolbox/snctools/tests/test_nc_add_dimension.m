function test_nc_add_dimension ( ncfile )
% TEST_NC_ADD_DIMENSION
%
% Relies upon nc_getdiminfo, nc_add_dimension.
%
% Test 1:  no inputs
% test 2:  too many inputs
% test 3:  first input not a netcdf file
% test 4:  2nd input not character
% test 5:  3rd input not numeric
% test 6:  3rd input is negative
% Test 7:  Add a normal length dimension.
% Test 8:  Add an unlimited dimension.
% test 9:  named dimension already exists (should fail)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_add_dimension.m 2306 2007-08-31 15:24:55Z johnevans007 $
% $LastChangedDate: 2007-08-31 11:24:55 -0400 (Fri, 31 Aug 2007) $
% $LastChangedRevision: 2306 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 'NC_ADD_DIMENSION:  starting test suite...\n' );

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
	nc_add_dimension;
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return





function test_002 ( ncfile )
% test 2:  too many inputs
testid = 'Test 2';
create_empty_file ( ncfile );
try
	nc_add_dimension ( ncfile, 'x', 10 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return









function test_003 ( ncfile )

% test 3:  first input not a netcdf file
create_empty_file ( ncfile );
fprintf ( 1, '    Skipping test 003.\n' );
%dimid = nc_add_dimension ( '/dev/null', 'x', 3 );
%fprintf ( 1, '%s:  passed.\n', testid );
return











function test_004 ( ncfile )

% test 4:  2nd input not char
create_empty_file ( ncfile );
try
	nc_add_dimension ( ncfile, 3, 3 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end

return










function test_005 ( ncfile )

% test 5:  3rd input not numeric
create_empty_file ( ncfile );
try
	dimid = nc_add_dimension ( ncfile, 't', 't' );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return









function test_006 ( ncfile )

% test 6:  3rd input is negative
create_empty_file ( ncfile );
try
	nc_add_dimension ( ncfile, 't', -1 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end
return







function test_007 ( ncfile )

% test 7:  add a normal dimension
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 5 );

%
% Now check that the new dimension are there.
d = nc_getdiminfo ( ncfile, 't' );
if ( ~strcmp(d.Name,'t') )
	error ( '%s:  nc_add_dimension failed on fixed dimension add name', mfilename  );
end
if ( d.Length ~= 5 )
	error ( '%s:  nc_add_dimension failed on fixed dimension add length', mfilename  );
end
if ( d.Unlimited ~= 0  )
	error ( '%s:  nc_add_dimension incorrectly classified the dimension', mfilename  );
end

return











function test_008 ( ncfile )
% test 8:  add an unlimited dimension
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );

%
% Now check that the new dimension are there.
d = nc_getdiminfo ( ncfile, 't' );
if ( ~strcmp(d.Name,'t') )
	error ( '%s:  nc_add_dimension failed on fixed dimension add name', mfilename  );
end
if ( d.Length ~= 0 )
	error ( '%s:  nc_add_dimension failed on fixed dimension add length', mfilename  );
end
if ( d.Unlimited ~= 1  )
	error ( '%s:  nc_add_dimension incorrectly classified the dimension', mfilename  );
end

return











function test_009 ( ncfile )

% test 9:  try to add a dimension that is already there
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
try
	nc_add_dimension ( ncfile, 't', 0 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
catch
	;
end
nc_add_dimension ( ncfile, 'x', 1 );
try
	nc_add_dimension ( ncfile, 'x', 1 );
	error ( '%s:  succeeded when it should have failed.\n', mfilename );
end

return





