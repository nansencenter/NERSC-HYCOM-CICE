function test_nc_datatype_string()
% TEST_NC_DATATYPE_STRING:
%
% Bad input argument tests.
% Test 1:  no inputs
% Test 2:  two inputs
% test 3:  input not numeric
% test 4:  input is outside of 0-6
%
% These tests should succeed
% test 5:  input is 0 ==> 'NC_NAT'
% test 6:  input is 1 ==> 'NC_BYTE'
% test 7:  input is 2 ==> 'NC_CHAR'
% test 8:  input is 3 ==> 'NC_SHORT'
% test 9:  input is 4 ==> 'NC_INT'
% test 10:  input is 5 ==> 'NC_FLOAT'
% test 11:  input is 6 ==> 'NC_DOUBLE'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_datatype_string.m 1808 2006-09-15 16:11:51Z johnevans007 $
% $LastChangedDate: 2006-09-15 12:11:51 -0400 (Fri, 15 Sep 2006) $
% $LastChangedRevision: 1808 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_DATATYPE_STRING: starting test suite ... \n' );

test_001;
test_002;
test_003;
test_004;
test_005;
test_006;
test_007;
test_008;
test_009;
test_010;
test_011;

return






function test_001 ( ncfile )

try
	dt = nc_datatype_string;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
catch
	;
end
return









function test_002 ( ncfile )

try
	dt = nc_datatype_string ( 0, 1 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
catch
	;
end
return










function test_003 ( ncfile )

% test 3:  input not numeric
try
	dt = nc_datatype_string ( 'a' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
catch
	;
end
return









function test_004 ( ncfile )


try
	dt = nc_datatype_string ( -1 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
catch
	;
end
return









function test_005 ( ncfile )
%
% These tests should succeed
% test 5:  input is 0 ==> 'NC_NAT'
dt = nc_datatype_string ( 0 );
if ~strcmp(dt,'NC_NAT')
	msg = sprintf ( '%s:  %s: failed to convert 0.\n', mfilename );
	error ( msg );
end
return










function test_006 ( ncfile )

dt = nc_datatype_string ( 1 );
if ~strcmp(dt,'NC_BYTE')
	msg = sprintf ( '%s:  %s: failed to convert 1.\n', mfilename );
	error ( msg );
end
return










function test_007 (ncfile )

dt = nc_datatype_string ( 2 );
if ~strcmp(dt,'NC_CHAR')
	msg = sprintf ( '%s:  %s: failed to convert 2.\n', mfilename );
	error ( msg );
end
return











function test_008 ( ncfile )

dt = nc_datatype_string ( 3 );
if ~strcmp(dt,'NC_SHORT')
	msg = sprintf ( '%s:  failed to convert 3.\n', mfilename );
	error ( msg );
end
return










function test_009 ( ncfile )

dt = nc_datatype_string ( 4 );
if ~strcmp(dt,'NC_INT')
	msg = sprintf ( '%s:  failed to convert 4.\n', mfilename );
	error ( msg );
end
return










function test_010 ( ncfile )

dt = nc_datatype_string ( 5 );
if ~strcmp(dt,'NC_FLOAT')
	msg = sprintf ( '%s:  failed to convert 5.\n', mfilename );
	error ( msg );
end
return










function test_011 ( ncfile )

dt = nc_datatype_string ( 6 );
if ~strcmp(dt,'NC_DOUBLE')
	msg = sprintf ( '%s: failed to convert 6.\n', mfilename );
	error ( msg );
end
return












