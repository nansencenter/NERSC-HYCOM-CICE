function test_nc_cat_a ( )
% TEST_NC_CAT_A:  tests the m-file nc_cat_a
%
% Test 01:  wrong number of input arguments.
% Test 02:  concatenating empty files
% Test 03:  First file has 5, 2nd has 10, 3rd has 15 time values.  They
%     do not overlap.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id$
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_CAT_A: starting test suite ... \n' );

test_01;
test_02;
test_03;

return;









%
% Create a netCDF file .  The general format is
%
% netcdf test_01 {
% dimensions:
%     time = UNLIMITED ; // (0 currently)
%     x = 10 ;
%     y = 20 ;
% variables:
%     double time(time) ;
%     double var1(time, y) ;
%     double var2(time, y, x) ;
% }
%
function create_default_test_file ( ncfile, abscissa_var );
nc_create_empty ( ncfile );

nc_add_dimension ( ncfile, 'time', 0 );
nc_add_dimension ( ncfile, 'x', 10 );
nc_add_dimension ( ncfile, 'y', 20 );

v.Name = 'time';
v.Dimension = { 'time' };
nc_addvar ( ncfile, v );

v.Name = 'var1';
v.Dimension = { 'time', 'y' };
nc_addvar ( ncfile, v );

v.Name = 'var2';
v.Dimension = { 'time', 'y', 'x' };
nc_addvar ( ncfile, v );

return








function test_01 ( );

ncfile1 = 'test_01.nc';
ncfile2 = 'test_02.nc';
ncfile3 = 'test_03.nc';
abscissa_var = 'time';

create_default_test_file ( ncfile1, abscissa_var );
create_default_test_file ( ncfile2, abscissa_var );
create_default_test_file ( ncfile3, abscissa_var );

try
	nc_cat_a;
	error ( '%s:  Succeeded when it should have failed.\n', mfilename );
end

return










function test_02 (  )

ncfile1 = 'test_01.nc';
ncfile2 = 'test_02.nc';
ncfile3 = 'test_03.nc';
abscissa_var = 'time';

create_default_test_file ( ncfile1, abscissa_var );
create_default_test_file ( ncfile2, abscissa_var );
create_default_test_file ( ncfile3, abscissa_var );

ncfiles{1} = 'test_01.nc';
ncfiles{2} = 'test_02.nc';
ncfiles{3} = 'test_03.nc';
output_ncfile = 'test_out.nc';
nc_cat_a ( ncfiles, output_ncfile, 'time' );
return










function test_03 (  )

ncfile1 = 'test_01.nc';
ncfile2 = 'test_02.nc';
ncfile3 = 'test_03.nc';
abscissa_var = 'time';

create_default_test_file ( ncfile1, abscissa_var );
vardata.time = [1:5]';
vardata.var1 = ones(5,20);
vardata.var2 = 2*ones(5,20,10);
nc_addnewrecs ( ncfile1, vardata, 'time' );



create_default_test_file ( ncfile2, abscissa_var );
vardata.time = [1:10]' + 10;
vardata.var1 = ones(10,20);
vardata.var2 = 2*ones(10,20,10);
nc_addnewrecs ( ncfile2, vardata, 'time' );



create_default_test_file ( ncfile3, abscissa_var );
vardata.time = [1:15]' + 20;
vardata.var1 = ones(15,20);
vardata.var2 = 2*ones(15,20,10);
nc_addnewrecs ( ncfile3, vardata, 'time' );


ncfiles = { ncfile1, ncfile2, ncfile3 };
output_ncfile = 'test_out.nc';
nc_cat_a ( ncfiles, output_ncfile, 'time' );

%
% There should be 30 records.
nt = nc_varsize ( output_ncfile, 'time' );
if nt ~= 30
	error ( '%s:  wrong number of final records.\n', mfilename );
end
return










