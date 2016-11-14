function test_snctools()
% TEST_SNCTOOLS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_snctools.m 2375 2007-10-23 12:08:08Z johnevans007 $
% $LastChangedDate: 2007-10-23 08:08:08 -0400 (Tue, 23 Oct 2007) $
% $LastChangedRevision: 2375 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'Pausing for ten seconds so you can read this...                    \n');
fprintf ( 1, '\n' );
fprintf ( 1, 'If you wish to test SNCTOOLS and OPeNDAP network access, you should\n' );
fprintf ( 1, 'set the SNCTOOLS TEST_REMOTE preference, i.e.                      \n' );
fprintf ( 1, '                                                                   \n' );
fprintf ( 1, '    setpref ( ''SNCTOOLS'', ''TEST_REMOTE'', true );               \n' );
fprintf ( 1, '                                                                   \n' );
fprintf ( 1, 'By default TEST_REMOTE is set to false since the remote servers    \n' );
fprintf ( 1, 'involved in the test may or may not be available.                  \n' );
fprintf ( 1, '                                                                   \n' );
pause(10);


mver = version('-release')
switch mver
	case { '2008a', '2007b', '2007a', '2006a', '2006b', '14', '13' }
		warning('off', 'SNCTOOLS:nc_archive_buffer:deprecatedMessage' );
		warning('off', 'SNCTOOLS:nc_datatype_string:deprecatedMessage' );
		warning('off', 'SNCTOOLS:nc_diff:deprecatedMessage' );
		warning('off', 'SNCTOOLS:nc_getall:deprecatedMessage' );
		warning('off', 'SNCTOOLS:snc2mat:deprecatedMessage' );
	case '12'
		error ( 'SNCTOOLS will not run under this version of matlab.' );
	otherwise
		warning ( 'This version of matlab has not been tested' );
end


%
% Save any old settings.
old_use_java_setting = getpref('SNCTOOLS','USE_JAVA',false);

setpref('SNCTOOLS','USE_JAVA',false);
run_all_tests;


% C indexing, java
setpref('SNCTOOLS','USE_JAVA',true);
run_all_tests;


fprintf ( 1, '\nAll SNCTOOLS tests succeeded.\n' );

fprintf ( 1, '\n' );
answer = input ( 'Do you wish to remove all test NetCDF and *.mat files that were created? [y/n]\n', 's' );
if strcmp ( lower(answer), 'y' )
    delete ( '*.nc' );
    delete ( '*.mat' );
end

fprintf ( 1, 'We''re done.\n' );




function run_all_tests()

test_nc_add_dimension    ( 'test.nc' );
test_nc_addhist          ( 'test.nc' );
test_nc_addvar           ( 'test.nc' );
test_nc_attget_attput    ( 'test.nc' );
test_nc_create_empty     ( 'test.nc' );
test_nc_datatype_string;
test_nc_iscoordvar       ( 'test.nc' );
test_nc_isunlimitedvar   ( 'test.nc' );
test_nc_isvar            ( 'test.nc' );
test_nc_varget_varput    ( 'test.nc' );
test_nc_varrename        ( 'test.nc' );
test_nc_varsize          ( 'test.nc' );

test_nc_dump             ( 'test.nc' );
test_nc_addnewrecs       ( 'test.nc' );
test_nc_diff             ( 'test1.nc', 'test2.nc' );
test_nc_add_recs         ( 'test.nc' );
test_nc_archive_buffer   ( 'test.nc' );
test_nc_cat_a;
test_nc_getall           ( 'test.nc' );
test_nc_getbuffer        ( 'test.nc' );
test_snc2mat             ( 'test.nc' );

test_nc_getdiminfo       ( 'test.nc' );
test_nc_getlast          ( 'test.nc' );
test_nc_getvarinfo       ( 'test.nc' );
test_nc_info             ( 'test.nc' );

return

