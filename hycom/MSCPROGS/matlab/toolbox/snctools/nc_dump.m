function nc_dump(file_name, varargin )
% NC_DUMP:  a Matlab counterpart to the NetCDF utility 'ncdump'.
%     NC_DUMP(NCFILE) prints metadata about the netCDF file NCFILE.  
%     NC_DUMP(NCFILE,VARNAME) prints metadata about just the one netCDF variable
%     named VARNAME.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_dump.m 2338 2007-09-21 18:41:58Z johnevans007 $
% $LastChangedDate: 2007-09-21 14:41:58 -0400 (Fri, 21 Sep 2007) $
% $LastChangedRevision: 2338 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


snc_nargchk(1,2,nargin);
snc_nargoutchk(0,0,nargout);


if nargin == 2
	do_restricted_variable = true;
	restricted_variable = varargin{1};
else
	do_restricted_variable = false;	
	restricted_variable = [];
end

metadata = nc_info ( file_name );


%
% print out name of file
fprintf ( 1, 'netcdf %s { \n\n', metadata.Filename );

dump_dimension_metadata ( metadata );

dump_variable_metadata ( metadata, restricted_variable );
	
if ( do_restricted_variable == false )
	dump_global_attributes ( metadata );
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump_dimension_metadata ( metadata )

if isfield ( metadata, 'Dimension' )
	num_dims = length(metadata.Dimension);
else
	num_dims = 0;
end

fprintf ( 1, 'dimensions:\n' );
for j = 1:num_dims
	if metadata.Dimension(j).Unlimited
		fprintf( 1, '\t%s = UNLIMITED ; (%i currently)\n', ...
		         deblank(metadata.Dimension(j).Name), metadata.Dimension(j).Length );
	else
		fprintf ( '\t%s = %i ;\n', metadata.Dimension(j).Name, metadata.Dimension(j).Length );
	end
end
fprintf('\n\n');

return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump_variable_metadata ( metadata, restricted_variable )

if isfield ( metadata, 'Dataset' )
	num_vars = length(metadata.Dataset);
else
	num_vars = 0;
end

fprintf ( 'variables:\n' );
for j = 1:num_vars

	if ~isempty(restricted_variable)
		if ~strcmp ( restricted_variable, metadata.Dataset(j).Name )
			continue
		end
	end

	dump_single_variable ( metadata.Dataset(j) );

end
fprintf ( '\n\n' );
return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump_single_variable ( var_metadata )

switch ( var_metadata.Nctype )
case 1
	fprintf ( 1, '\tbyte ' );
case 2
	fprintf ( 1, '\tchar ' );
case 3
	fprintf ( 1, '\tshort ' );
case 4
	fprintf ( 1, '\tlong ' );
case 5
	fprintf ( 1, '\tfloat ' );
case 6
	fprintf ( 1, '\tdouble ' );
end

fprintf ( 1, '%s', var_metadata.Name );

if isempty(var_metadata.Dimension) 
	fprintf ( 1, '([]), ' );
else
	fprintf ( 1, '(%s', var_metadata.Dimension{1} );
	for k = 2:length(var_metadata.Size)
		fprintf ( 1, ',%s', var_metadata.Dimension{k} );
	end
	fprintf ( 1, '), ');
end


if isempty(var_metadata.Dimension)
	fprintf ( 1, 'shape = [1]\n' );
else
	fprintf ( 1, 'shape = [%d', var_metadata.Size(1)  );
	for k = 2:length(var_metadata.Size)
		fprintf ( 1, ' %d', var_metadata.Size(k)  );
	end
	fprintf ( 1, ']\n');
end


%
% Now do all attributes for each variable.
num_atts = length(var_metadata.Attribute);
for k = 1:num_atts
	dump_single_attribute ( var_metadata.Attribute(k), var_metadata.Name );
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump_single_attribute ( attribute, varname )

switch ( attribute.Nctype )
case 0
	att_val = '';
	att_type = 'NC_NAT';
case 1
	att_val = sprintf ('%d ', fix(attribute.Value) );
	att_type = 'x';
case 2
	att_val = sprintf ('"%s" ', attribute.Value );
	att_type = '';
case 3
	att_val = sprintf ('%i ', attribute.Value );
	att_type = 's';
case 4
	att_val = sprintf ('%i ', attribute.Value );
	att_type = 'd';
case 5
	att_val = sprintf ('%f ', attribute.Value );
	att_type = 'f';
case 6
	att_val = sprintf ('%g ', attribute.Value );
	att_type = '';
end
if nargin == 1
	fprintf( 1, '\t\t:%s = %s%s\n', ...
         attribute.Name, att_val, att_type);
else
	fprintf( 1, '\t\t%s:%s = %s%s\n', ...
         varname, attribute.Name, att_val, att_type);
end

return
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump_global_attributes ( metadata )


if isfield ( metadata, 'Attribute' )
	num_atts = length(metadata.Attribute);
else
	num_atts = 0;
end

if num_atts > 0
	fprintf ( 1, '//global attributes:\n' );
end

for k = 1:num_atts
	dump_single_attribute ( metadata.Attribute(k) );
end


fprintf ( 1, '}\n' );

return
