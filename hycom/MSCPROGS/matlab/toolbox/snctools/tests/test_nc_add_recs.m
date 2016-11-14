function test_nc_add_recs ( ncfile )
% TEST_NC_ADD_RECS
%
% Relies upon nc_getvarino, nc_addvar
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
%   10.  Do two successive writes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_add_recs.m 1808 2006-09-15 16:11:51Z johnevans007 $
% $LastChangedDate: 2006-09-15 12:11:51 -0400 (Fri, 15 Sep 2006) $
% $LastChangedRevision: 1808 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ADD_RECS:  starting test suite...\n' );

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

return





function create_ncfile ( ncfile );

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
varstruct.Name = 'test_var3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };

nc_addvar ( ncfile, varstruct );

return





function test_001 ( ncfile )

%
% Try no inputs
try
	nc_add_recs;
	error ( [ mfilename ':' ' nc_add_recs succeeded on no inputs, should have failed'] );
end
return








function test_002 ( ncfile )
%
% Try one inputs
test_id = 'Test 2';
try
	nc_add_recs ( ncfile );
	error ( '%s:  nc_add_recs succeeded on one input, should have failed', mfilename );
end
return








function test_003 ( ncfile )

%
% Try with 2nd input that isn't a structure.
test_id = 'Test 3';
try
	nc_add_recs ( ncfile, [] );
	error ( '%s:  nc_add_recs succeeded on one input, should have failed', mfilename );
end
return











function test_004 ( ncfile )

%
% Try with 2nd input that is an empty structure.
test_id = 'Test 4';
try
	nc_add_recs ( ncfile, struct([]) );
	error ( '%s:  nc_add_recs succeeded on empty structure, should have failed', mfilename );
end
return











function test_005 ( ncfile )

%
% Try a structure with bad names
test_id = 'Test 5';
input_data.a = [3 4];
input_data.b = [5 6];
try
	nc_add_recs ( ncfile, input_data );
	error ( '%s:  nc_add_recs succeeded on a structure with bad names, should have failed', mfilename );
end
return







function test_006 ( ncfile )

%
% Try good data with a bad unlimited dimension name
test_id = 'Test 6';
input_data.test_var = [3 4]';
input_data.test_var2 = [5 6]';
try
	nc_add_recs ( ncfile, input_data, 'bad_unlim_dim_name' );
	error ( '%s:  nc_add_recs succeeded with a badly named unlimited dimension, should have failed', mfilename );
end
return






function test_007 ( ncfile )



%
% Try a good test.
test_id = 'Test 7';
before = nc_getvarinfo ( ncfile, 'test_var2' );


input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';

nc_add_recs ( ncfile, input_buffer );

after = nc_getvarinfo ( ncfile, 'test_var2' );
if ( (after.Size - before.Size) ~= 3 )
	error ( '%s:  nc_add_recs failed to add the right number of records.', mfilename );
end


return











function test_008 ( ncfile )

%
% Try writing to a fixed size variable
test_id = 'Test 8';


input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';
input_buffer.test_var3 = [3 4 5]';

try
	nc_add_recs ( ncfile, input_buffer );
	error ( '%s:  nc_add_recs succeeded on writing to a fixed size variable, should have failed.', mfilename );
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
	nc_add_recs ( ncfile, input_buffer );
	error ( '%s:  nc_add_recs passed when writing to a file with no unlimited dimension', mfilename );
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
varstruct.Name = 'test_var3';
varstruct.Nctype = 'double';
varstruct.Dimension = { 'x' };

nc_addvar ( ncfile, varstruct );


before = nc_getvarinfo ( ncfile, 'test_var2' );
clear input_buffer;
input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';
nc_add_recs ( ncfile, input_buffer );
nc_add_recs ( ncfile, input_buffer );

after = nc_getvarinfo ( ncfile, 'test_var2' );
if ( (after.Size - before.Size) ~= 6 )
	error ( '%s:  nc_add_recs failed to add the right number of records.', mfilename );
end
return











