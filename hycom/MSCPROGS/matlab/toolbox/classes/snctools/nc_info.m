function fileinfo = nc_info ( ncfile )
% NC_INFO:  information about a NetCDF 2 or 3 file
%
% USAGE:  fileinfo = nc_info ( ncfile );
%
% PARAMETERS:
% Input:
%    ncfile:  
%        a string that specifies the name of the NetCDF file
% Output:
%    fileinfo:
%        A structure whose fields contain information about the contants
%        of the NetCDF file.  The set of fields return in "fileinfo" are:
%
%        Filename:  
%            a string containing the name of the file.
%        Dimension:  
%            an array of structures describing the NetCDF dimensions.
%        Dataset:  
%            an array of structures describing the NetCDF datasets.
%        Attributes:  
%            An array of structures These correspond to the global attributes.
%
%
%        Each "Dimension" element contains the following fields:
%       
%        Name:
%            a string containing the name of the dimension.
%        Length:
%            a scalar value, the size of this dimension
%        Unlimited:
%            Set to 1 if the dimension is the record dimension, set to
%            0 otherwise.
%
%
%        Each "Dataset" element contains the following structures.
%
%        Name:  
%            a string containing the name of the variable.
%        Nctype:  
%            a string specifying the NetCDF datatype of this variable.
%        Dimensions:  
%            a cell array with the names of the dimensions upon which
%            this variable depends.
%        Unlimited:  
%            Flag, either 1 if the variable has an unlimited dimension
%            or 0 if not.
%        Rank:  
%            Array that describes the size of each dimension upon which 
%            this dataset depends.
%        DataAttributes:  
%            Same as "Attributes" above, but here they are the variable 
%            attributes.
%                         
%        Each "Attribute" or "DataAttribute" element contains the following 
%        fields.
%
%        Name:  
%            a string containing the name of the attribute.
%        Nctype:  
%            a string specifying the NetCDF datatype of this attribute.
%        Attnum:  
%            a scalar specifying the attribute id
%        Value: 
%            either a string or a double precision value corresponding to
%            the value of the attribute
%
%
% The "Dataset" elements are not populated with the actual data values.
%
% This routine purposefully mimics that of Mathwork's hdfinfo.
%
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_info.m 2340 2007-09-21 18:46:07Z johnevans007 $
% $LastChangedDate: 2007-09-21 14:46:07 -0400 (Fri, 21 Sep 2007) $
% $LastChangedRevision: 2340 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



snc_nargchk(1,1,nargin);
snc_nargoutchk(1,1,nargout);


%
% Do we use java instead of mexnc?
use_java = getpref ( 'SNCTOOLS', 'USE_JAVA', 0 );
if use_java 
	fileinfo = nc_info_java(ncfile);
else
	fileinfo = nc_info_mex(ncfile);
end

return











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileinfo = nc_info_mex ( ncfile )


fileinfo.Filename = ncfile;

[ncid, status]=mexnc('open', ncfile, nc_nowrite_mode );
if status ~= 0
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_INFO:MEXNC:OPEN', ncerr );
end



[ndims, nvars, ngatts, record_dimension, status] = mexnc('INQ', ncid);
if status ~= 0
    ncerr = mexnc('strerror', status);
    mexnc('close',ncid);
    snc_error ( 'SNCTOOLS:NC_INFO:MEXNC:INQ', ncerr );
end


%
% Get the dimensions
if ndims == 0
	Dimension = struct ( [] );
else
	if ndims > 0
		Dimension(1)=nc_getdiminfo ( ncid, 0 );
	end
	Dimension = repmat ( Dimension, ndims,1 );
	for dimid = 1:ndims-1
		Dimension(dimid+1)=nc_getdiminfo ( ncid, dimid );
	end
end



%
% Get the global attributes.
if ngatts == 0
	fileinfo.Attribute = struct([]);
else
	if ngatts > 0
		Attribute(1) = nc_get_attribute_struct ( ncid, nc_global, 0 );
	end
	Attribute = repmat ( Attribute, ngatts, 1 );
	for attnum = 1:ngatts-1
		Attribute(attnum+1) = nc_get_attribute_struct ( ncid, nc_global, attnum );
	end
	fileinfo.Attribute = Attribute;
end





%
% Get the variable information.
if nvars == 0
	Dataset = struct([]);
else
	if ( nvars > 0 )
		Dataset(1) = nc_getvarinfo ( ncid, 0 );
	end
	Dataset = repmat ( Dataset, nvars, 1 );
	for varid=1:nvars-1
		Dataset(varid+1) = nc_getvarinfo ( ncid, varid );
	end
end

fileinfo.Dimension = Dimension;
fileinfo.Dataset = Dataset;


mexnc('close',ncid);


return









