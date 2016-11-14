function new_data = nc_archive_buffer ( input_buffer, ncfile, record_variable )
% NC_ARCHIVE_BUFFER:  Tacks on new data from simple matlab structure to an unlimited-dimension netcdf file
% 
% This function is deprecated.  Please use nc_addnewrecs.m instead.  All 
% this function really does is call nc_addnewrecs anyway.
%
% If a field is present in the structure, but not in the netcdf file, then that
% field is ignored.  Only data that is more recent than the last record currently
% in the NetCDF file is written.  This routine is written with time series data in mind.
%
% USAGE:  new_data = nc_archive_buffer ( buffer, ncfile, record_variable, unlimited_dimension )
% 
% PARAMETERS:
%   Input:
%      buffer:  
%          structure of (hopefully) time series data.  
%      ncfile:  
%          netcdf file that we write information to
%      record_variable:
%          hopefully this is something like "time".  In ROMS, it is "ocean_time".
%   Output:
%      new_data:  
%          Matlab structure of data corresponding in structure to "buffer", but
%          consisting only of those records which were actually written to file.
%
% In case of an error, an exception is thrown.
%
%  
%   The dimensions of the data should match that of the target netcdf file.  For example, 
%   suppose an ncdump of the
%   NetCDF file looks something like
%
%       netcdf a_netcdf_file {
%       dimensions:
%       	lat = 1 ;
%       	lon = 2 ;
%       	depth = 2 ; 
%       	time = UNLIMITED ; // (500 currently)
%       variables:
%       	double time(time) ;
%       	float var1(time, depth) ;
%       	float var2(time, depth, lat, lon) ;
%       	float var3(time, depth, lat) ;
%       
%       // global attributes:
%       }
% 
%   The "input_buffer" should look something like the following:
%
%       >> input_buffer
%
%       input_buffer =
%
%           time: [3x1 double]
%           var1: [3x2 double]
%           var2: [4-D double]
%           var3: [3x2 double]
%
%   The reason for the possible size discrepency here is that matlab will
%   ignore trailing singleton dimensions (but not interior ones, such as
%   that in var2.
%
%   If a variable is not present, then the corresponding NetCDF variable will
%   populate with the appropriate _FillValue.
%          
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_archive_buffer.m 1870 2006-10-09 12:08:53Z johnevans007 $
% $LastChangedDate: 2006-10-09 08:08:53 -0400 (Mon, 09 Oct 2006) $
% $LastChangedRevision: 1870 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wid = sprintf ( 'SNCTOOLS:%s:deprecatedMessage', lower(mfilename) );
msg = sprintf( '%s is deprecated and may be removed in a future version of SNCTOOLS.', upper(mfilename) );
warning ( wid, msg );

new_data = [];

snc_nargchk(2,3,nargin);

switch nargin
case 2
	new_data = nc_addnewrecs ( ncfile, input_buffer );
case { 3, 4 }
	new_data = nc_addnewrecs ( ncfile, input_buffer, record_variable );
end


return;



