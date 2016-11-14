function test_nc_archive_buffer ( ncfile )
% TEST_NC_ADDNEWRECS
%
% Relies upon nc_addvar, nc_getvarinfo
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_archive_buffer.m 1808 2006-09-15 16:11:51Z johnevans007 $
% $LastChangedDate: 2006-09-15 12:11:51 -0400 (Fri, 15 Sep 2006) $
% $LastChangedRevision: 1808 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 'NC_ARCHIVE_BUFFER:  starting test suite...\n' );


create_test_ncfile ( ncfile );
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
test_011 ( ncfile );


return










function create_test_ncfile ( ncfile );

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
	nb = nc_archive_buffer;
	error ( [ mfilename ':' ' nc_archive_buffer succeeded on no inputs, should have failed'] );
end
return










function test_002 ( ncfile )

%
% Try one inputs
try
	nb = nc_archive_buffer ( [] );
	error ( [ mfilename ': succeeded on one input, should have failed'] );
end
return










function test_003 ( ncfile )

%
% Try with 2nd input that isn't a structure.
try
	nb = nc_archive_buffer ( [], ncfile );
	error ( [ mfilename ': succeeded, should have failed'] );
end
return









function test_004 ( ncfile )

%
% Try with 2nd input that is an empty structure.
try
	nb = nc_archive_buffer ( struct([]), ncfile );
	error ( [ mfilename ': succeeded, should have failed'] );
end
return










function test_005 ( ncfile )

%
% Try a structure with bad names
test_id = 'Test 5';
input_data.a = [3 4];
input_data.b = [5 6];
try
	nb = nc_archive_buffer ( input_data, ncfile );
	error ( [ mfilename ': succeeded, should have failed'] );
end
return





function test_006 ( ncfile )

%
% Try good data with a bad record variable name
input_data.test_var = [3 4]';
input_data.test_var2 = [5 6]';
try
	nd = nc_archive_buffer ( input_data, ncfile, 'bad_time' );
	error ( [ mfilename ': succeeded, should have failed'] );
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

new_data = nc_archive_buffer ( input_buffer, ncfile, 'time' );

after = nc_getvarinfo ( ncfile, 'test_var2' );
if ( (after.Size - before.Size) ~= 3 )
	error ( 'nc_archive_buffer failed to add the right number of records.' );
end
return










function test_008 ( ncfile )

%
% Try writing to a fixed size variable


input_buffer.test_var = single([3 4 5]');
input_buffer.test_var2 = [3 4 5]';
input_buffer.test_var3 = [3 4 5]';

try
	nc_archive_buffer ( input_buffer, ncfile );
	error ( 'nc_archive_buffer succeeded on writing to a fixed size variable, should have failed.' );
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
	nb = nc_archive_buffer ( input_buffer,  ncfile );
	error ( 'nc_archive_buffer passed when writing to a file with no unlimited dimension' );
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


new_data = nc_archive_buffer ( input_buffer, ncfile );
input_buffer.time = [4 5 6]';
new_data = nc_archive_buffer ( input_buffer, ncfile );

after = nc_getvarinfo ( ncfile, 'time' );
if ( (after.Size - before.Size) ~= 6 )
	error ( 'nc_archive_buffer failed to add the right number of records.' );
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


new_data = nc_archive_buffer ( input_buffer, ncfile );
input_buffer.time = [3 4 5]';
new_data = nc_archive_buffer ( input_buffer, ncfile );

after = nc_getvarinfo ( ncfile, 'time' );
if ( (after.Size - before.Size) ~= 5 )
	error ( 'nc_archive_buffer failed to add the right number of records.' );
end


return
