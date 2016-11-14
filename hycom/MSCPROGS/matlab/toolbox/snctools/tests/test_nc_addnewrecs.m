function test_nc_addnewrecs ( ncfile )
% TEST_NC_ADDNEWRECS
%
% Relies on nc_addvar, nc_getvarinfo
%
% Test run include
%    1.  No inputs, should fail.
%    2.  One inputs, should fail.
%    3.  Two inputs, 2nd is not a structure, should fail.
%    4.  Two inputs, 2nd is an empty structure, should fail.
%    5.  Two inputs, 2nd is a structure with bad variable names, should fail.
%    6.  Three inputs, 3rd is non existant unlimited dimension.
%    7.  Two inputs, write to two variables, should succeed.
%    8.  Two inputs, write to two variables, one of them not unlimited, should fail.
%    9.  Try to write to a file with no unlimited dimension.
%   10.  Do two successive writes.  Should succeed.
%   11.  Do two successive writes, but on the 2nd write let the coordinate
%        variable overlap with the previous write.  Should still succeed,
%        but fewer datums will be written out.
%   12.  Do two successive writes, but with the same data.  Should 
%        return an empty buffer, but not fail
% Test 13:  Add a single record.  This is a corner case.
% Test 14:  Add a single record, trailing singleton dimensions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_addnewrecs.m 2264 2007-08-08 19:27:49Z johnevans007 $
% $LastChangedDate: 2007-08-08 15:27:49 -0400 (Wed, 08 Aug 2007) $
% $LastChangedRevision: 2264 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ADDNEWRECS:  starting test suite...\n' );

create_ncfile ( ncfile )

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








function create_ncfile ( ncfile )

%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end


%
% Create a fixed dimension.  
len_x = 4;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 )
	error ( mexnc('strerror',status) );
end

%
% Create two singleton dimensions.
[ydimid, status] = mexnc ( 'def_dim', ncid, 'y', 1 );
if ( status ~= 0 )
	error ( mexnc('strerror',status) );
end

[zdimid, status] = mexnc ( 'def_dim', ncid, 'z', 1 );
if ( status ~= 0 )
	error ( mexnc('strerror',status) );
end



len_t = 0;
[ydimid, status] = mexnc ( 'def_dim', ncid, 'time', 0 );
if ( status ~= 0 )
	error ( mexnc('strerror',status) );
end





%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end

%
% Add a variable along the time dimension
varstruct.Name = 'test_var';
varstruct.Nctype = 'float';
varstruct.Dimension = { 'time' };
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(1).Value = 'This is a test';
varstruct.Attribute(2).Name = 'short_val';
varstruct.Attribute(2).Value = int16(5);

nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 'test_var2';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };

nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 'trailing_singleton';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time', 'z', 'y' };

nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 'time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };

nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 'test_var3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };

nc_addvar ( ncfile, varstruct );

return










function test_001 ( ncfile )

%
% Try no inputs
try
	nc_addnewrecs;
	error ( [ mfilename ':' ' nc_addnewrecs succeeded on no inputs, should have failed'] );
end
return







function test_002 ( ncfile )

%
% Try one inputs
try
	nc_addnewrecs ( ncfile );
	error ( '%s:  nc_addnewrecs succeeded on one input, should have failed', mfilename );
end
return











function test_003 ( ncfile )

%
% Try with 2nd input that isn't a structure.
try
	nc_addnewrecs ( ncfile, [] );
	error ( '%s:  nc_addnewrecs succeeded on one input, should have failed', mfilename );
end
return












function test_004 ( ncfile )

%
% Try with 2nd input that is an empty structure.
try
	nc_addnewrecs ( ncfile, struct([]) );
	error ( '%s:  nc_addnewrecs succeeded on empty structure, should have failed', mfilename );
end
return











function test_005 ( ncfile )

%
% Try a structure with bad names
input_data.a = [3 4];
input_data.b = [5 6];
try
	nb = nc_addnewrecs ( ncfile, input_data );
	error ( '%s:  nc_addnewrecs succeeded on a structure with bad names, should have failed', mfilename );
end
return











function test_006 ( ncfile )

%
% Try good data with a bad record variable name
input_data.test_var = [3 4]';
input_data.test_var2 = [5 6]';
try
	nb = nc_addnewrecs ( ncfile, input_data, 'bad_time' );
	error ( '%s:  nc_addnewrecs succeeded with a badly named record variable, should have failed', mfilename );
end
return









function test_007 ( ncfile )

%
% Try a good test.
before = nc_getvarinfo ( ncfile, 'test_var2' );


clear input_buffer;
input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';
input_buffer.time = [1 2 3]';

new_data = nc_addnewrecs ( ncfile, input_buffer );

after = nc_getvarinfo ( ncfile, 'test_var2' );
if ( (after.Size - before.Size) ~= 3 )
	error ( '%s:  nc_addnewrecs failed to add the right number of records.', mfilename );
end
return











function test_008 ( ncfile )

%
% Try writing to a fixed size variable


input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';
input_buffer.test_var3 = [3 4 5]';

try
	nc_addnewrecs ( ncfile, input_buffer );
	error ( '%s:  nc_addnewrecs succeeded on writing to a fixed size variable, should have failed.', mfilename );
end

return












function test_009 ( ncfile )

%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end


%
% Create a fixed dimension.  
len_x = 4;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim x, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end

clear varstruct;
varstruct.Name = 'test_var3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };

nc_addvar ( ncfile, varstruct );


input_buffer.time = [1 2 3]';
try
	nb = nc_addnewrecs ( ncfile, input_buffer );
	error ( '%s:  nc_addnewrecs passed when writing to a file with no unlimited dimension', mfilename );
end
return











function test_010 ( ncfile )

%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end


%
% Create a fixed dimension.  
len_x = 4;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim x, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


len_t = 0;
[ydimid, status] = mexnc ( 'def_dim', ncid, 'time', 0 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim time, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end





%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end

clear varstruct;
varstruct.Name = 'time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };
nc_addvar ( ncfile, varstruct );




before = nc_getvarinfo ( ncfile, 'time' );

clear input_buffer;
input_buffer.time = [1 2 3]';


new_data = nc_addnewrecs ( ncfile, input_buffer );
input_buffer.time = [4 5 6]';
new_data = nc_addnewrecs ( ncfile, input_buffer );

after = nc_getvarinfo ( ncfile, 'time' );
if ( (after.Size - before.Size) ~= 6 )
	error ( '%s:  nc_addnewrecs failed to add the right number of records.', mfilename );
end

return








function test_011 ( ncfile )

%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end


%
% Create a fixed dimension.  
len_x = 4;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim x, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


len_t = 0;
[ydimid, status] = mexnc ( 'def_dim', ncid, 'time', 0 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim time, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end





%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end

clear varstruct;
varstruct.Name = 'time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };
nc_addvar ( ncfile, varstruct );




before = nc_getvarinfo ( ncfile, 'time' );

clear input_buffer;
input_buffer.time = [1 2 3]';


new_data = nc_addnewrecs ( ncfile, input_buffer );
input_buffer.time = [3 4 5]';
new_data = nc_addnewrecs ( ncfile, input_buffer );

after = nc_getvarinfo ( ncfile, 'time' );
if ( (after.Size - before.Size) ~= 5 )
	error ( '%s:  nc_addnewrecs failed to add the right number of records.', mfilename );
end
return














function test_012 ( ncfile )

%
% baseline case
create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 'ocean_time', 0 );

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
nb = nc_addnewrecs ( ncfile, b, 'ocean_time' );
if ( ~isempty(nb) )
	msg = sprintf ( '%s:  nc_addnewrecs failed on %s.\n', mfilename, ncfile );
	error ( msg );
end
v = nc_getvarinfo ( ncfile, 't1' );
if ( v.Size ~= 10 )
	error ( '%s:  expected var length was not 10.\n', mfilename );
end

return








function test_013 ( ncfile )

create_empty_file ( ncfile );

nc_add_dimension ( ncfile, 'time', 0 );
nc_add_dimension ( ncfile, 'x', 10 );
nc_add_dimension ( ncfile, 'y', 10 );

clear varstruct;
varstruct.Name = 'time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't1';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };
nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 't2';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time', 'y', 'x' };
nc_addvar ( ncfile, varstruct );


clear varstruct;
varstruct.Name = 't3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time', 'y' };
nc_addvar ( ncfile, varstruct );


b.time = [0];
b.t1 = 0;
b.t2 = zeros(10,10);
b.t3 = zeros(1,10);

nc_addnewrecs ( ncfile, b, 'time' );

clear b
b.time = [1];
b.t1 = 1;
b.t2 = ones(10,10);
b.t3 = ones(1,10);
nc_addnewrecs ( ncfile, b, 'time' );


%
% Now read them back.  
b = nc_getbuffer ( ncfile, 0, 2 );
if length(b.time) ~= 2
	msg = sprintf ( '%s:  length of time variable was %d and not the expected 2\n', mfilename, length(b.time) );
	error ( msg );
end
if (b.time(1) ~= 0) & (b.time(2) ~= 1)
	msg = sprintf ( '%s:  values of time variable are wrong\n', mfilename );
	error ( msg );
end


return









function test_014 ( ncfile )

create_empty_file ( ncfile );

nc_add_dimension ( ncfile, 'time', 0 );
nc_add_dimension ( ncfile, 'x', 1 );
nc_add_dimension ( ncfile, 'y', 1 );

clear varstruct;
varstruct.Name = 'time';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time' };
nc_addvar ( ncfile, varstruct );

clear varstruct;
varstruct.Name = 't1';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'time', 'y', 'x' };
nc_addvar ( ncfile, varstruct );


b.time = [0];
b.t1 = 0;

nc_addnewrecs ( ncfile, b, 'time' );

clear b
b.time = [1];
b.t1 = 1;
nc_addnewrecs ( ncfile, b, 'time' );


%
% Now read them back.  
t1 = nc_varget ( ncfile, 't1' );
if (t1(1) ~= 0) & (t1(2) ~= 1)
	msg = sprintf ( '%s:  values are wrong\n', mfilename );
	error ( msg );
end


return
