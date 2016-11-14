function test_nc_create_empty ( ncfile )
%
% Test 01:  Supply no arguments.
% Test 02:  No mode given
% Test 03:  64-bit mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id$
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ( 1, 'NC_CREATE_EMPTY: starting test suite ... \n' );

test_01 ( ncfile );
test_02 ( ncfile );
test_03 ( ncfile );

return





function test_01 ( ncfile )

try
	nc_create_empty;
	msg = sprintf ( '%s:  succeeded when it should have failed\n', mfilename );
	error ( msg );
end

return








function test_02 ( ncfile )

nc_create_empty ( ncfile );
md = nc_info ( ncfile );

if length(md.Dataset) ~= 0
	msg = sprintf ( '%s:  number of variables was not zero\n', mfilename );
	error ( msg );
end

if length(md.Attribute) ~= 0
	msg = sprintf ( '%s:  number of global attributes was not zero\n', mfilename );
	error ( msg );
end

if length(md.Dimension) ~= 0
	msg = sprintf ( '%s:  number of dimensions was not zero\n', mfilename );
	error ( msg );
end

return







function test_03 ( ncfile )

mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );

nc_create_empty ( ncfile, mode );
md = nc_info ( ncfile );

if length(md.Dataset) ~= 0
	msg = sprintf ( '%s:  number of variables was not zero\n', mfilename );
	error ( msg );
end

if length(md.Attribute) ~= 0
	msg = sprintf ( '%s:  number of global attributes was not zero\n', mfilename );
	error ( msg );
end

if length(md.Dimension) ~= 0
	msg = sprintf ( '%s:  number of dimensions was not zero\n', mfilename );
	error ( msg );
end

return


