function Dataset = nc_getvarinfo_java ( arg1, arg2 )
%
% This function handles the java case.


switch ( class(arg1) )
case { 'ucar.nc2.NetcdfFile', 'ucar.nc2.Group', 'ucar.nc2.dods.DODSNetcdfFile' }
	Dataset = snc_java_varid_info ( arg2 );
case 'char'
	Dataset = get_varinfo_closed ( arg1, arg2 );
end

return





%===============================================================================
function Dataset = get_varinfo_closed ( ncfile, varname )

import ucar.nc2.dods.*     % import opendap reader classes
import ucar.nc2.*          % have to import this (NetcdfFile) as well for local reads
                           
if snc_is_url ( ncfile )
	try  
		jncid = DODSNetcdfFile(ncfile);
	catch
		try
			jncid = NetcdfFile.open(ncfile);
		catch
			snc_error ( 'SNCTOOLS:nc_getvarinfo_java:urlOpen', ...
			        'Could not open netCDF URL.' );
		end
	end
else
	try
		jncid = NetcdfFile.open(ncfile);
	catch
		snc_error ( 'SNCTOOLS:nc_getvarinfo_java:localFileOpen', ...
		        'Could not open netCDF file.' );
	end
end

jvarid = jncid.findVariable(varname);
if isempty(jvarid)
	close(jncid);
	msg = sprintf ('Could not locate variable %s', varname );
	snc_error ( 'SNCTOOLS:NC_GETVARINFO:badVariableName', msg );
end



Dataset = snc_java_varid_info ( jvarid );

close ( jncid );

return


%===============================================================================
function Dataset = get_varinfo_open ( jncid, jvarid )




%
% All the details are hidden here because we need the exact same
% functionality in nc_info.


return
