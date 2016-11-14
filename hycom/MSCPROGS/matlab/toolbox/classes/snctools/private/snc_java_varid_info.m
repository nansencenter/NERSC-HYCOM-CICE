function DataSet = snc_java_varid_info ( jvarid )
% SNC_JAVA_VARID_INFO:  returns metadata structure for a netcdf variable
%
% This function is private to SNCTOOLS.  It is called by nc_info and
% nc_getvarinfo, and uses the java API.
%
% USAGE:   DataSet = snc_java_varid_info ( jvarid );
% 
% PARAMETERS:
% Input:
%     jvarid:  
%         of type ucar.nc2.dods.DODSVariable
% Output:
%     DataSet:
%         array of metadata structures.  The fields are
%         
%         Name
%         Nctype
%         Unlimited
%         Dimension
%         Attribute
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id$
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DataSet.Name = char ( jvarid.getName() );

%
% Get the datatype, store as an integer
datatype = char(jvarid.getDataType().toString());
switch ( datatype )
case 'double'
	DataSet.Nctype = nc_double;
case 'float'
	DataSet.Nctype = nc_float;
case 'int'
	DataSet.Nctype = nc_int;
case 'short'
	DataSet.Nctype = nc_short;

%
% So apparently, DODSNetcdfFile returns 'String', while
% NetcdfFile returns 'char'???
case { 'String', 'char' }
	DataSet.Nctype = nc_char;
case 'byte'
	DataSet.Nctype = nc_byte;
otherwise 
	msg = sprintf ( '%s:  unhandled datatype ''%s''\n', mfilename, datatype );
	error ( msg );
end




%
% determine if it is unlimited or not
DataSet.Unlimited = double ( jvarid.isUnlimited() );


%
% Retrieve the dimensions
Dimension = {};
dims = jvarid.getDimensions();
nvdims = dims.size();
for j = 1:nvdims
	theDim = jvarid.getDimension(j-1);
	Dimension{j} = char ( theDim.getName() );
end
DataSet.Dimension = Dimension;



%
% Get the size of the variable
if nvdims == 0
	DataSet.Size = 1;
else
	Size = double ( jvarid.getShape() );
	DataSet.Size = Size';
end

%
% Get the list of attributes.
%
Attribute = [];
j_att_list = jvarid.getAttributes();
DataSet.Attribute = snc_java_bundle_atts ( j_att_list );


return

