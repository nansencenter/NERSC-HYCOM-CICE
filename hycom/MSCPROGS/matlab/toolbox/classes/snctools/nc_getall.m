function data = nc_getall ( ncfile )
% NCDATA = NC_GETALL(NCFILE) reads the entire contents of the netCDF file 
% NCFILE into the structure NCDATA.  This function is deprecated and may be
% dropped in a future release of SNCTOOLS.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_getall.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wid = sprintf ( 'SNCTOOLS:%s:deprecatedMessage', lower(mfilename) );
msg = sprintf( '%s is deprecated and may be removed in a future version of SNCTOOLS.', upper(mfilename) );
warning ( wid, msg );

data = [];


% Show usage if too few arguments.
%
if nargin~=1 
	error ( 'must have one input argument.\n' );
end


%
% Open netCDF file
%
[cdfid,status ]=mexnc('open',ncfile,'NOWRITE');
if status ~= 0
	error ( mexnc('strerror', status) );
end



[nvars, status] = mexnc('INQ_NVARS', cdfid);
if status < 0
	mexnc('close',cdfid);
	error ( mexnc('strerror', status) );
end
[ngatts, status] = mexnc('INQ_NATTS', cdfid);
if status < 0
	mexnc('close',cdfid);
	error ( mexnc('strerror', status) );
end

for varid=0:nvars-1

	varstruct = [];

	[varname, datatype, ndims, dims, natts, status] = mexnc('INQ_VAR', cdfid, varid);
	if status < 0 
		mexnc('close',cdfid);
		error ( mexnc('strerror', status) );
	end


	%
	% If ndims is zero, then it must be a singleton variable.  Don't bother trying
	% to retrieve the data, there is none.
	if ( ndims == 0 )
		varstruct.data = [];
	else
		values = nc_varget(ncfile, varname);
		varstruct.data = values;
	end



	%
	% get all the attributes
	for attnum = 0:natts-1

		[attname, status] = mexnc('inq_attname', cdfid, varid, attnum);
		if status < 0 
			mexnc('close',cdfid);
			error ( mexnc('strerror', status) );
		end

		try
			attval = nc_attget(ncfile, varname, attname);
		catch
			mexnc('close',cdfid);
			msg2 = sprintf ( 'nc_attget failed, ''%s''.\n', lasterr );
			error ( msg2 );
		end
		
		%
		% Matlab structures don't like the leading '_'
		if strcmp(attname,'_FillValue' )
			attname = 'FillValue';
		end


		sanitized_attname = matlab_sanitize_attname ( attname );


		%
		% this puts the attribute into the variable structure
		varstruct.(sanitized_attname) = attval;


	end


	%
	% Add this variable to the entire file structure
	data.(varname) = varstruct;

end


%
% Do the same for the global attributes
%
% get all the attributes
varname = 'global';
global_atts = [];
for attnum = 0:ngatts-1

	[attname, status] = mexnc('inq_attname', cdfid, nc_global, attnum);
	if status < 0 
		mexnc('close',cdfid);
		error ( mexnc('strerror',status) );
	end

	try
		attval = nc_attget(ncfile, nc_global, attname);
	catch
		mexnc('close',cdfid);
		msg = sprintf ( 'nc_attget failed, ''%s''.\n', lasterr );
		error ( msg );
	end
	
	sanitized_attname = matlab_sanitize_attname ( attname );


	%
	% this puts the attribute into the variable structure
	global_atts.(sanitized_attname) = attval;



end

if ~isempty ( global_atts )
	data.global_atts = global_atts;
end

mexnc('close',cdfid);

if isempty(data)
	data = struct([]);
end

return



function sanitized_attname = matlab_sanitize_attname ( attname )
	%
	% could the attribute name  be interpreted as a number?
	% If so, must fix this.
	% An attribute name of, say, '0' is not permissible in matlab
	if ~isnan(str2double(attname))
		sanitized_attname = ['SNC_' attname];
	elseif ( attname(1) == '_' )
		sanitized_attname = ['SNC_' attname];
	else
		sanitized_attname = attname;

		%
		% Does the attribute have non-letters in the leading 
		% position?  Convert these to underscores.  
		if ( ~isletter(attname(1)) );
			sanitized_attname(1) = '_';
		end
	end
return


