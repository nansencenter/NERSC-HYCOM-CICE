function test_nc_getbuffer ( ncfile )
% TEST_NC_GETBUFFER
%
% Relies upon nc_addvar, nc_addnewrecs, nc_add_dimension
%
% test 1:  no input arguments, should fail
% test 2:  2 inputs, 2nd is not a cell array, should fail
% test 3:  3 inputs, 2nd and 3rd are not numbers, should fail
% test 4:  4 inputs, 2nd is not a cell array, should fail
% test 5:  4 inputs, 3rd and 4th are not numbers, should fail
% test 6:  1 input, 1st is not a file, should fail.
% test 7:  5 inputs, should fail
% test 8:  1 input, an empty netcdf with no variables, should fail
%          because no record variable was found
%
% test 9:  1 input, 3 record variables. Should succeed.
% test 10:  2 inputs, same netcdf file as 9.  Restrict output to two
%           of the variables.  Should succeed.
% test 11:  3 inputs, same netcdf file as 9.  Restrict output to given
%           start:start+count range, which is given as valid.
% test 12:  3 inputs, same netcdf file as 9.  Restrict output to given
%           start:start+count range.  Start is negative number.  Result 
%           should be the last few "count" records.
% test 13:  3 inputs, same netcdf file as 9.  Restrict output to given
%           start:start+count range.  count is negative number.  Result 
%           should be everything from start to "end - count"
% test 14:  4 inputs.  Otherwise the same as test 11.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_getbuffer.m 1807 2006-09-14 11:04:27Z johnevans007 $
% $LastChangedDate: 2006-09-14 07:04:27 -0400 (Thu, 14 Sep 2006) $
% $LastChangedRevision: 1807 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_GETBUFFER:  starting test suite...\n' );
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
try
	nb = nc_getbuffer;
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_002 ( ncfile )
try
	nb = nc_getbuffer ( ncfile, 1 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_003 ( ncfile )
try
	nb = nc_getbuffer ( ncfile, 'a', 'b' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_004 ( ncfile )
try
	nb = nc_getbuffer ( ncfile, 1, 1, 2 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_005 ( ncfile )
try
	nb = nc_getbuffer ( ncfile, cell(1), 'a', 'b' );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_006 ( ncfile )
try
	nb = nc_getbuffer ( 5 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_007 ( ncfile )
try
	nb = nc_getbuffer ( ncfile, cell(1), 3, 4, 5 );
	msg = sprintf ( '%s:   succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return





function test_008 ( ncfile )
%
% empty file case
create_empty_file ( ncfile );
try
	nb = nc_getbuffer ( ncfile );
	msg = sprintf ( '%s:  : succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end
return




function test_009 ( ncfile )
%
% baseline case
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

nb = nc_getbuffer ( ncfile );

%
% should have 4 fields
f = fieldnames(nb);
n = length(f);
if n ~= 4
	msg = sprintf ( '%s:  : output buffer did not have 4 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:4
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 10 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end
return





function test_010 ( ncfile )
%
% baseline case
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

nb = nc_getbuffer ( ncfile );

%
% restrict variables

nb = nc_getbuffer ( ncfile, {'t1', 't2'} );

%
% should have 2 fields
f = fieldnames(nb);
n = length(f);
if n ~= 2
	msg = sprintf ( '%s:  : output buffer did not have 2 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:2
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 10 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end
return






function test_011 ( ncfile )


%
% baseline case
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

nb = nc_getbuffer ( ncfile );

nb = nc_getbuffer ( ncfile, 5, 3 );

%
% should have 4 fields
f = fieldnames(nb);
n = length(f);
if n ~= 4
	msg = sprintf ( '%s:  : output buffer did not have 4 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:n
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 3 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end

%
% t1 should be [5 6 7]
if any ( nb.t1 - [5 6 7]' )
	msg = sprintf ( '%s:  : t1 was not what we thought it should be.\n', mfilename  );
	error ( msg );
end
return







function test_012 ( ncfile )

%
% baseline case
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

%
% start is negative

nb = nc_getbuffer ( ncfile, -1, 3 );

%
% should have 4 fields
f = fieldnames(nb);
n = length(f);
if n ~= 4
	msg = sprintf ( '%s:  : output buffer did not have 4 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:n
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 3 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end

%
% t1 should be [7 8 9]
if any ( nb.t1 - [7 8 9]' )
	msg = sprintf ( '%s:  : t1 was not what we thought it should be.\n', mfilename  );
	error ( msg );
end
return







function test_013 ( ncfile )

%
% baseline case
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

nb = nc_getbuffer ( ncfile, 5, -1 );
%
% should have 4 fields
f = fieldnames(nb);
n = length(f);
if n ~= 4
	msg = sprintf ( '%s:  : output buffer did not have 4 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:n
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 5 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end

%
% t1 should be [5 6 7 8 9]
if any ( nb.t1 - [5 6 7 8 9]' )
	msg = sprintf ( '%s:  : t1 was not what we thought it should be.\n', mfilename  );
	error ( msg );
end
return







function test_014 ( ncfile )


%
% baseline case
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

nb = nc_getbuffer ( ncfile, {'t1', 't2' }, 5, -1 );

%
% should have 2 fields
f = fieldnames(nb);
n = length(f);
if n ~= 2
	msg = sprintf ( '%s:  : output buffer did not have 2 fields.\n', mfilename  );
	error ( msg );
end
for j = 1:n
	fname = f{j};
	d = getfield ( nb, fname );
	if ( length(d) ~= 5 )
		msg = sprintf ( '%s:  : length of field %s in the output buffer was not 10.\n', mfilename  );
		error ( msg );
	end
end

%
% t1 should be [5 6 7 8 9]
if any ( nb.t1 - [5 6 7 8 9]' )
	msg = sprintf ( '%s:  : t1 was not what we thought it should be.\n', mfilename  );
	error ( msg );
end
return




