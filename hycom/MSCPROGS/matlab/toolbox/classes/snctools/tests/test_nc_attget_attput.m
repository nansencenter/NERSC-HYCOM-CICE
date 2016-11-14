function test_nc_attget_attput ( ncfile )
% TEST_NC_ATTGET_ATTPUT
%
% Tests run include
%
% 1.  write/retrieve a double attribute
% 2.  write/retrieve a float attribute
% 3.  write/retrieve a int attribute
% 4.  write/retrieve a short int attribute
% 5.  write/retrieve a uint8 attribute
% 6.  write/retrieve a int8 attribute
% 7.  write/retrieve a text attribute
% 9.  write/retrieve a global attribute, using -1 as the variable name.
% 10.  write/retrieve a global attribute, using nc_global as the variable name.
% 11.  write/retrieve a global attribute, using 'GLOBAL' as the variable name.
% 12.  Try to retrieve a non existing attribute, should fail
% 21.  write/retrieve a new double attribute
% 22.  write/retrieve a new float attribute
% 23.  write/retrieve a new int attribute
% 24.  write/retrieve a new short int attribute
% 25.  write/retrieve a new uint8 attribute
% 26.  write/retrieve a new int8 attribute
% 27.  write/retrieve a new text attribute

% 401:  try to retrieve an attribute from a non dods url

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_attget_attput.m 2368 2007-10-11 13:48:43Z johnevans007 $
% $LastChangedDate: 2007-10-11 09:48:43 -0400 (Thu, 11 Oct 2007) $
% $LastChangedRevision: 2368 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ATTGET, NC_ATTPUT:  starting test suite...\n' );


create_test_ncfile ( ncfile )

test_01 ( ncfile );
test_02 ( ncfile );
test_03 ( ncfile );
test_04 ( ncfile );
test_05 ( ncfile );
test_06 ( ncfile );
test_07 ( ncfile );
test_08 ( ncfile );
test_09 ( ncfile );
test_10 ( ncfile );
test_11 ( ncfile );
test_12 ( ncfile );

test_21 ( ncfile );
test_22 ( ncfile );
test_23 ( ncfile );
test_24 ( ncfile );
test_25 ( ncfile );
test_26 ( ncfile );
test_27 ( ncfile );



return




function test_401 (ncfile)

if ~usejava('jvm' )
	return
end
if ~ ( getpref ( 'SNCTOOLS', 'TEST_REMOTE', false ) )
	return;
end
url = 'http://rocky.umeoce.maine.edu/GoMPOM/cdfs/gomoos.20070723.cdf';
fprintf ( 1, 'Testing remote URL access %s...\n', url );
w = nc_attget ( url, 'w', 'valid_range' );
if ~strcmp(class(w),'single')
	error ( 'Class of retrieve attribute was not single' );
end
if (abs(double(w(1)) - 0.5) < eps)
	error ( 'valid max did not match' );
end
if (abs(double(w(2)) - 0.5) < eps)
	error ( 'valid max did not match' );
end
return


function test_03 ( ncfile )

attvalue = nc_attget ( ncfile, 'x_db', 'test_int_att' );
if ( ~strcmp(class(attvalue), 'int32' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int32.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= int32(3) )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename);
	error ( msg );
end

return











function test_04 ( ncfile )




attvalue = nc_attget ( ncfile, 'x_db', 'test_short_att' );
if ( ~strcmp(class(attvalue), 'int16' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int16.\n', mfilename );
	error ( msg );
end
if ( length(attvalue) ~= 2 )
	msg = sprintf ( '%s:  retrieved attribute length differs from what was written.\n', mfilename );
	error ( msg );
end
if ( any(double(attvalue) - [5 7])  )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename  );
	error ( msg );
end

return






function test_05 ( ncfile )

attvalue = nc_attget ( ncfile, 'x_db', 'test_uchar_att' );
if ( ~strcmp(class(attvalue), 'int8' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int8.\n', mfilename );
	error ( msg );
end
if ( uint8(attvalue) ~= uint8(100) )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return




function test_06 ( ncfile )

attvalue = nc_attget ( ncfile, 'x_db', 'test_schar_att' );
if ( ~strcmp(class(attvalue), 'int8' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int8.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= int8(-100) )
	msg = sprintf ( '%s:  %s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return






function test_07 ( ncfile )

attvalue = nc_attget ( ncfile, 'x_db', 'test_text_att' );
if ( ~strcmp(class(attvalue), 'char' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not char.\n', mfilename );
	error ( msg );
end

if ( ~strcmp(attvalue,'abcdefghijklmnopqrstuvwxyz') )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return





function test_08 ( ncfile )

warning ( 'off', 'SNCTOOLS:nc_attget:java:doNotUseGlobalString' );

attvalue = nc_attget ( ncfile, '', 'test_double_att' );
if ( ~strcmp(class(attvalue), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= 3.14159 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

warning ( 'on', 'SNCTOOLS:nc_attget:java:doNotUseGlobalString' );

return





function test_09 ( ncfile )

attvalue = nc_attget ( ncfile, -1, 'test_double_att' );
if ( ~strcmp(class(attvalue), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= 3.14159 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return





function test_10 ( ncfile )

attvalue = nc_attget ( ncfile, nc_global, 'test_double_att' );
if ( ~strcmp(class(attvalue), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= 3.14159 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename  );
	error ( msg );
end

return 






function test_11 ( ncfile )

warning ( 'off', 'SNCTOOLS:nc_attget:doNotUseGlobalString' );
warning ( 'off', 'SNCTOOLS:nc_attget:java:doNotUseGlobalString' );

attvalue = nc_attget ( ncfile, 'GLOBAL', 'test_double_att' );
if ( ~strcmp(class(attvalue), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= 3.14159 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename  );
	error ( msg );
end

warning ( 'on', 'SNCTOOLS:nc_attget:java:doNotUseGlobalString' );
warning ( 'on', 'SNCTOOLS:nc_attget:doNotUseGlobalString' );

return





function test_12 ( ncfile )

try
	attvalue = nc_attget ( ncfile, 'z_double', 'test_double_att' );
	msg = sprintf ( '%s:  %s:  nc_attget succeeded when it should have failed.\n', mfilename  );
	error ( msg );
end

return















function create_test_ncfile ( ncfile )
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


len_y = 6;
[ydimid, status] = mexnc ( 'def_dim', ncid, 'y', len_y );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_dim'' failed on dim y, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


[xvarid, status] = mexnc ( 'def_var', ncid, 'x_db', 'NC_DOUBLE', 2, [ydimid xdimid] );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_var'' failed on var x_short, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


[zvarid, status] = mexnc ( 'def_var', ncid, 'z_double', 'NC_DOUBLE', 2, [ydimid xdimid] );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''def_var'' failed on var x_short, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end

%
% Define attributes for all datatypes for x_db, but not z_double
% The short int attribute will have length 2
status = mexnc ( 'put_att_double', ncid, xvarid, 'test_double_att', nc_double, 1, 3.14159 );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_float', ncid, xvarid, 'test_float_att', nc_float, 1, single(3.14159) );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_int', ncid, xvarid, 'test_int_att', nc_int, 1, int32(3) );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_short', ncid, xvarid, 'test_short_att', nc_short, 2, int16([5 7]) );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_uchar', ncid, xvarid, 'test_uchar_att', nc_byte, 1, uint8(100) );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_schar', ncid, xvarid, 'test_schar_att', nc_byte, 1, int8(-100) );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_text', ncid, xvarid, 'test_text_att', nc_char, 26, 'abcdefghijklmnopqrstuvwxyz' );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end
status = mexnc ( 'put_att_double', ncid, nc_global, 'test_double_att', nc_double, 1, 3.14159 );
if ( status < 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  put_att_double failed, ''%s''\n', mfilename, ncerr_msg );
	error ( msg );
end




[status] = mexnc ( 'end_def', ncid );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  ''end_def'' failed, file %s, error message '' %s ''\n', mfilename, ncfile, ncerr_msg );
	error ( msg );
end


%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	error ( 'CLOSE failed' );
end


return





function test_01 ( ncfile )

attvalue = nc_attget ( ncfile, 'x_db', 'test_double_att' );
if ( ~strcmp(class(attvalue), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end
if ( attvalue ~= 3.14159 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return




function test_02 ( ncfile )

create_test_ncfile ( ncfile );

attvalue = nc_attget ( ncfile, 'x_db', 'test_float_att' );
if ( ~strcmp(class(attvalue), 'single' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not single.\n', mfilename );
	error ( msg );
end
if ( abs(double(attvalue) - 3.14159) > 1e-6 )
	msg = sprintf ( '%s:  retrieved attribute differs from what was written.\n', mfilename );
	error ( msg );
end

return



function test_21 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', 0 );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'double' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not double.\n', mfilename );
	error ( msg );
end

if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end

return




function test_22 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', single(0) );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'single' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not single.\n', mfilename );
	error ( msg );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


function test_23 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', int32(0) );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'int32' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int32.\n', mfilename );
	error ( msg );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


function test_24 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', int16(0) );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'int16' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int16.\n', mfilename );
	error ( msg );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


function test_25 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', int8(0) );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'int8' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not int8.\n', mfilename );
	error ( msg );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


function test_26 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', uint8(0) );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'int8' ) )
	msg = sprintf ( 'class of retrieved attribute was %s and not int8.\n', class(x) );
	error ( msg );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


function test_27 ( ncfile )

nc_attput ( ncfile, 'x_db', 'new_att', '0' );
x = nc_attget ( ncfile, 'x_db', 'new_att' );

if ( ~strcmp(class(x), 'char' ) )
	msg = sprintf ( '%s:  class of retrieved attribute was not char.\n', mfilename );
	error ( msg );
end
if (x ~= '0' )
	error ( 'retrieved attribute was not same as written value' );
end


return





