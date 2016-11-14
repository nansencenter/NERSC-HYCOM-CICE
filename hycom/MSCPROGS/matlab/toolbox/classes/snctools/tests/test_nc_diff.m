function test_nc_diff ( ncfile_1, ncfile_2 )
% TEST_NC_DIFF
%
% Relies upon nc_addnewrecs, nc_addvar
%
% Test run include
% 1.  Create two empty identical files.
% 2.  Add different data to the test_var variable in both files. 
% 3.  Test files with unequal dimensions.
% 4.  Test files with variables with unequal rank.
% 5.  Test files and attributes.
% 6.  Test files and attributes with a different number of attributes.  
%     Should fail.
% 7.  Test files and attributes with a different number of attributes.  
%     Reverse order.  Should fail.
% 8.  Test files and attributes, one attribute has different class.
% 9.  Test files and attributes, one attribute has different length
% 10.  Test files and attributes, one attribute has different value
% 11.  Test files and global attributes, one attribute has different value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_diff.m 1808 2006-09-15 16:11:51Z johnevans007 $
% $LastChangedDate: 2006-09-15 12:11:51 -0400 (Fri, 15 Sep 2006) $
% $LastChangedRevision: 1808 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 'NC_DIFF:  starting test suite...\n' );

test_001 ( ncfile_1, ncfile_2 );
test_002 ( ncfile_1, ncfile_2 );
test_003 ( ncfile_1, ncfile_2 );
test_004 ( ncfile_1, ncfile_2 );
test_005 ( ncfile_1, ncfile_2 );
test_006 ( ncfile_1, ncfile_2 );
test_007 ( ncfile_1, ncfile_2 );
test_008 ( ncfile_1, ncfile_2 );
test_009 ( ncfile_1, ncfile_2 );
test_010 ( ncfile_1, ncfile_2 );
test_011 ( ncfile_1, ncfile_2 );
return





% 
% Test the empty files
function test_001 ( ncfile_1, ncfile_2 )

create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );

status = nc_diff ( ncfile_1, ncfile_2 );
if  status < 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return






function test_002 ( ncfile_1, ncfile_2 )

create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );

write_buffer.time = [1:5]';
write_buffer.test_var = rand(5,1);
nd = nc_addnewrecs ( ncfile_1, write_buffer, 'time' );

write_buffer.test_var = rand(5,1);
nd = nc_addnewrecs ( ncfile_2, write_buffer, 'time' );

status = nc_diff ( ncfile_1, ncfile_2 );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return




function test_003 ( ncfile_1, ncfile_2 )

create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );

write_buffer.time = [1:5]';
write_buffer.test_var = rand(5,1);
nd = nc_addnewrecs ( ncfile_1, write_buffer, 'time' );

write_buffer.test_var = rand(5,1);
nd = nc_addnewrecs ( ncfile_2, write_buffer, 'time' );

%
% add to the unlimited dimension in just one file
write_buffer.time = [6:7]';
write_buffer.test_var = rand(2,1);
nd = nc_addnewrecs ( ncfile_1, write_buffer, 'time' );

status = nc_diff ( ncfile_1, ncfile_2 );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return










function test_004 ( ncfile_1, ncfile_2 )

create_test_file ( ncfile_1 );
write_buffer.time = [1:5]';
write_buffer.test_var = rand(5,1);
nd = nc_addnewrecs ( ncfile_1, write_buffer, 'time' );

%
% add to the unlimited dimension in just one file
write_buffer.time = [6:7]';
write_buffer.test_var = rand(2,1);
nd = nc_addnewrecs ( ncfile_1, write_buffer, 'time' );
%
% Create a netcdf file with a variable with different rank
create_test_file ( ncfile_2, 'rank different' );
status = nc_diff ( ncfile_1, ncfile_2 );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end

return










function test_005 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and attributes.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
status = nc_diff ( ncfile_1, ncfile_2, '-attributes' );
if  status < 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return










function test_006 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and different number of variable attributes.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_1, 'test_var', 'test_attribute', 0 );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return









function test_007 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and different number of variable attributes.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_2, 'test_var', 'test_attribute', 0 );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return















function test_008 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and one attribute with different class.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_1, 'test_var', 'short_val', 0 );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return






function test_009 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and one attribute with different class.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_1, 'test_var', 'short_val', int16([5 3]) );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return








function test_010 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and one attribute with different value.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_1, 'test_var', 'short_val', int16([27]) );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end
return






function test_011 ( ncfile_1, ncfile_2 )

% 
% Test the empty files and one global attribute with different value.
create_test_file ( ncfile_1 );
create_test_file ( ncfile_2 );
nc_attput ( ncfile_1, nc_global, 'short_val', int16([27]) );
status = nc_diff ( ncfile_1, ncfile_2, '-attribute' );
if  status >= 0
	msg = sprintf ( '%s:   did not produce the expected result.\n', mfilename  );
	error ( msg );
end

return








function create_test_file ( ncfile, arg2 )

%
% ok, first create the first file
[ncid_1, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''create'' failed, error message '' %s ''\n', mfilename, ncerr_msg );
	error ( msg );
end


%
% Create a fixed dimension.  
len_x = 4;
[xdimid, status] = mexnc ( 'def_dim', ncid_1, 'x', len_x );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim x, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end

%
% Create a fixed dimension.  
len_y = 5;
[ydimid, status] = mexnc ( 'def_dim', ncid_1, 'y', len_y );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim y, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


len_t = 0;
[ydimid, status] = mexnc ( 'def_dim', ncid_1, 'time', 0 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim time, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end





%
% CLOSE
status = mexnc ( 'close', ncid_1 );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end

%
% Add a variable along the time dimension
varstruct.Name = 'test_var';
varstruct.Nctype = 'float';
if nargin > 1
	varstruct.Dimension = { 'time', 'y' };
else
	varstruct.Dimension = { 'time' };
end
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


