function Dataset = nc_getvarinfo ( arg1, arg2 )
% NC_GETVARINFO:  returns metadata about a specific NetCDF variable
%
% VINFO = NC_GETVARINFO(NCFILE,VARNAME) returns a metadata structure VINFO about
% the variable VARNAME in the netCDF file NCFILE.
%
% VINFO = NC_GETVARINFO(NCID,VARID) returns a metadata structure VINFO about
% the variable whose netCDF variable-id is VARID, and whose parent file-id is 
% NCID.  The netCDF file is assumed to be open, and in this case the file will
% not be closed upon completion.
%
% VINFO will have the following fields:
%
%    Name:  
%       a string containing the name of the variable.
%    Nctype:  
%       a string specifying the NetCDF datatype of this variable.
%    Unlimited:  
%       Flag, either 1 if the variable has an unlimited dimension or 0 if not.
%    Dimensions:  
%       a cell array with the names of the dimensions upon which this variable 
%       depends.
%    Attribute:  
%       An array of structures corresponding to the attributes defined for the 
%       specified variable.
%                         
%    Each "Attribute" element contains the following fields.
%
%       Name:  
%           a string containing the name of the attribute.
%       Nctype:  
%           a string specifying the NetCDF datatype of this attribute.
%       Attnum:  
%           a scalar specifying the attribute id
%       Value: 
%           either a string or a double precision value corresponding to the 
%           value of the attribute
%
% In case of an error, an exception is thrown.
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_getvarinfo.m 2353 2007-10-03 21:45:10Z johnevans007 $
% $LastChangedDate: 2007-10-03 17:45:10 -0400 (Wed, 03 Oct 2007) $
% $LastChangedRevision: 2353 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Show usage if too few arguments.
%
snc_nargchk(2,2,nargin);
snc_nargoutchk(1,1,nargout);




%
% Do we use java instead of mexnc?
use_java = getpref ( 'SNCTOOLS', 'USE_JAVA', false );
if use_java 
	Dataset = nc_getvarinfo_java(arg1,arg2);
else
	Dataset = nc_getvarinfo_mex(arg1,arg2);
end

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dataset = nc_getvarinfo_mex ( arg1, arg2 )

%
% If we are here, then we must have been given something local.
if ischar(arg1) && ischar(arg2)

	ncfile = arg1;
	varname = arg2;


	[ncid,status ]=mexnc('open',ncfile,nc_nowrite_mode);
	if status ~= 0
    	ncerr = mexnc('strerror', status);
	    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:OPEN', ncerr );
	end


	[varid, status] = mexnc('INQ_VARID', ncid, varname);
	if ( status ~= 0 )
    	ncerr = mexnc('strerror', status);
	    mexnc('close',ncid);
	    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ncerr );
	end

	
	Dataset = get_varinfo ( ncid,  varid );

	%
	% close whether or not we were successful.
	mexnc('close',ncid);


elseif isnumeric ( arg1 ) && isnumeric ( arg2 )

	ncid = arg1;
	varid = arg2;

	Dataset = get_varinfo ( ncid,  varid );

else
	snc_error ( 'SNCTOOLS:NC_GETVARINFO:badTypes', ...
	        'Must have either both character inputs, or both numeric.' );
end


return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dataset = get_varinfo ( ncid, varid )


[record_dimension, status] = mexnc ( 'INQ_UNLIMDIM', ncid );
if status ~= 0
   	ncerr = mexnc('strerror', status);
    mexnc('close',ncid);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_UNLIMDIM', ncerr );
end



[varname, datatype, ndims, dims, natts, status] = mexnc('INQ_VAR', ncid, varid);
if status ~= 0 
   	ncerr = mexnc('strerror', status);
    mexnc('close',ncid);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_VAR', ncerr );
end



Dataset.Name = varname;
Dataset.Nctype = datatype;

%
% Assume the current variable does not have an unlimited dimension until
% we know that it does.
Dataset.Unlimited = false;

if ndims == 0
	Dataset.Dimension = {};
	Dataset.Size = 1;
else

	for j = 1:ndims
	
		[dimname, dimlength, status] = mexnc('INQ_DIM', ncid, dims(j));
		if ( status ~= 0 )
   			ncerr = mexnc('strerror', status);
		    mexnc('close',ncid);
		    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_DIM', ncerr );
		end
	
		Dataset.Dimension{j} = dimname; 
		Dataset.Size(j) = dimlength;
	
		if dims(j) == record_dimension
			Dataset.Unlimited = true;
		end
	end
end

%
% get all the attributes
if natts == 0
	Dataset.Attribute = struct([]);
else
	for attnum = 0:natts-1
		Dataset.Attribute(attnum+1) = nc_get_attribute_struct ( ncid, varid, attnum );
	end
end







return











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dataset = nc_getvarinfo_java ( ncfile, varname )
%
% This function handles the java case.

import ucar.nc2.dods.*     % import opendap reader classes
import ucar.nc2.*          % have to import this (NetcdfFile) as well for local reads
                           


%
% Open netCDF file
%
if snc_is_url ( ncfile )
		jncid = DODSNetcdfFile(ncfile);
else
	jncid = NetcdfFile(ncfile);
end


jvarid = jncid.findVariable(varname);
if isempty(jvarid)
	close(jncid);
	msg = sprintf ('Could not locate variable %s', varname );
	error ( 'SNCTOOLS:NC_GETVARINFO:badVariableName', msg );
end


%
% All the details are hidden here because we need the exact same
% functionality in nc_info.
Dataset = snc_java_varid_info ( jvarid );


close ( jncid );

return
