function nc_varrename ( ncfile, old_variable_name, new_variable_name )
% NC_VARRENAME:  renames a NetCDF variable.
%
% NC_VARRENAME(NCFILE,OLD_VARNAME,NEW_VARNAME) renames a netCDF variable from
% OLD_VARNAME to NEW_VARNAME.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_varrename.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


snc_nargchk(3,3,nargin);
snc_nargoutchk(0,0,nargout);


[ncid,status ]=mexnc('OPEN',ncfile,nc_write_mode);
if status ~= 0
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:OPEN', ncerr );
end


status = mexnc('REDEF', ncid);
if status ~= 0
	mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:REDEF', ncerr );
end


[varid, status] = mexnc('INQ_VARID', ncid, old_variable_name);
if status ~= 0
	mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ncerr );
end


status = mexnc('RENAME_VAR', ncid, varid, new_variable_name);
if status ~= 0
	mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:RENAME_VAR', ncerr );
end


status = mexnc('ENDDEF', ncid);
if status ~= 0
	mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:ENDDEF', ncerr );
end


status = mexnc('close',ncid);
if status ~= 0
    ncerr = mexnc('strerror', status);
    snc_error ( 'SNCTOOLS:NC_VARGET:MEXNC:CLOSE', ncerr );
end

return;
