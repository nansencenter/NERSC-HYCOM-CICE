function test_nc_addhist ( ncfile )
% TEST_NC_ADDHIST
%
% Relies upon nc_attget, nc_add_dimension, nc_addvar
% Test 1:  no inputs
% test 2:  too many inputs
% test 3:  first input not a netcdf file
% test 4:  2nd input not character
% test 5:  3rd input not character
% Test 6:  Add history first time to global attributes
% Test 7:  Add history again
% Test 8:  Add history first time to variable
% Test 9:  Add history again
% Test 10:  Add history to a variable.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_nc_addhist.m 2311 2007-08-31 20:35:46Z johnevans007 $
% $LastChangedDate: 2007-08-31 16:35:46 -0400 (Fri, 31 Aug 2007) $
% $LastChangedRevision: 2311 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_ADDHIST:  starting test suite\n' );

test_01 ( ncfile );
test_02 ( ncfile );
test_03 ( ncfile );
test_04 ( ncfile );
%test_05 ( ncfile );
test_06 ( ncfile );
test_07 ( ncfile );
%test_08 ( ncfile );
%test_09 ( ncfile );
%test_10 ( ncfile );


return


function test_01 ( ncfile );
try
	nc_addhist;
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return



function test_02 ( ncfile );
try
	nc_addhist ( ncfile, 'x', 'blurb', 'blurb' );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return






function test_03 ( ncfile )

create_empty_file ( ncfile );
try
	nc_addhist ( 'asdfjsadjfsadlkjfsa;ljf;l', 'test' );
	msg = sprintf ( '%s:  %s succeeded when it should have failed.\n', mfilename, testid );
	error ( msg );
end
return





function test_04 ( ncfile )

create_empty_file ( ncfile );
try
	nc_addhist ( ncfile, 5 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return




function test_05 ( ncfile )

create_empty_file ( ncfile );
nc_add_dimension ( ncfile, 't', 0 );
clear varstruct;
varstruct.Name = 'T';
nc_addvar ( ncfile, varstruct );
try
	nc_addhist ( ncfile, 'T', 5 );
	msg = sprintf ( '%s:  succeeded when it should have failed.\n', mfilename );
	error ( msg );
end
return





function test_06 ( ncfile )

create_empty_file ( ncfile );
histblurb = 'blah';
nc_addhist ( ncfile, histblurb );

hista = nc_attget ( ncfile, nc_global, 'history' );
s = findstr(hista, histblurb );
if isempty(hista)
	msg = sprintf ( '%s:  history attribute did not contain first attribution.\n', mfilename );
	error ( msg );
end
return




function test_07 ( ncfile )

create_empty_file ( ncfile );
histblurb = 'blah a';
nc_addhist ( ncfile, histblurb );
histblurb2 = 'blah b';
nc_addhist ( ncfile, histblurb2 );
histatt = nc_attget ( ncfile, nc_global, 'history' );
s = findstr(histatt, histblurb2 );
if isempty(histatt)
	msg = sprintf ( '%s:  history attribute did not contain second attribution.\n', mfilename );
	error ( msg );
end
return




function test_08 ( ncfile )

create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'T';
nc_addvar ( ncfile, varstruct );
histblurb = 'blah';
nc_addhist ( ncfile, 'T', histblurb );
hista = nc_attget ( ncfile, 'T', 'history' );
s = findstr(hista, histblurb );
if isempty(hista)
	msg = sprintf ( '%s:  history attribute did not contain first attribution.\n', mfilename );
	error ( msg );
end
return






function test_09 ( ncfile )

create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'T';
nc_addvar ( ncfile, varstruct );
histblurb = 'blah a';
nc_addhist ( ncfile, 'T', histblurb );
histblurb2 = 'blah b';
nc_addhist ( ncfile, 'T', histblurb2 );
histatt = nc_attget ( ncfile, 'T', 'history' );
s = findstr(histatt, histblurb2 );
if isempty(histatt)
	msg = sprintf ( '%s:  history attribute did not contain second attribution.\n', mfilename );
	error ( msg );
end
return












function test_10 ( ncfile )

create_empty_file ( ncfile );
clear varstruct;
varstruct.Name = 'x';
nc_addvar ( ncfile, varstruct );
histblurb = 'blah a';
nc_addhist ( ncfile, 'x', histblurb );
histatt = nc_attget ( ncfile, 'x', 'history' );
if ~strcmp ( histatt(end-5:end), histblurb )
	msg = sprintf ( '%s:  history attribute did not match the input.\n', mfilename );
	error ( msg );
end
return













