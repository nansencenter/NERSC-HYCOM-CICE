function nc_addvar ( ncfile, varstruct )
% NC_ADDVAR:  adds a variable to a NetCDF file
%
% USAGE:  nc_addvar ( ncfile, varstruct );
%
% PARAMETERS:
% Input
%    ncfile:
%    varstruct:
%        This is a structure with four fields:
%
%        Name
%        Nctype
%        Dimension
%        Attribute
%
%      "Name" is just that, the name of the variable to be defined.
%
%      "Nctype" should be 
%          'double', 'float', 'int', 'short', or 'byte', or 'char'
%          'NC_DOUBLE', 'NC_FLOAT', 'NC_INT', 'NC_SHORT', 'NC_BYTE', 'NC_CHAR'
%
%      "Dimension" is a cell array of dimension names.
%
%      "Attribute" is also a structure array.  Each element has two
%      fields, "Name", and "Value".
%
% Output: 
%     None.  In case of an error, an exception is thrown.
%
% AUTHOR:
%    john.g.evans.ne@gmail.com
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_addvar.m 2340 2007-09-21 18:46:07Z johnevans007 $
% $LastChangedDate: 2007-09-21 14:46:07 -0400 (Fri, 21 Sep 2007) $
% $LastChangedRevision: 2340 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snc_nargchk(2,2,nargin);


if  ~ischar(ncfile) 
	snc_error ( 'SNCTOOLS:NC_ADDVAR:badInput', 'file argument must be character' );
end

if ( ~isstruct(varstruct) )
	snc_error ( 'SNCTOOLS:NC_ADDVAR:badInput', '2nd argument must be a structure' );
end


varstruct = validate_varstruct ( varstruct );



[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	msg = sprintf ( 'OPEN failed on %s, ''%s''', ncfile, ncerr );
	snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:OPEN', msg );
end

%
% determine the dimids of the named dimensions
num_dims = length(varstruct.Dimension);
dimids = zeros(num_dims,1);
for j = 1:num_dims
	[dimids(j), status] = mexnc ( 'dimid', ncid, varstruct.Dimension{j} );
	if ( status ~= 0 )
		mexnc ( 'close', ncid );
		ncerr = mexnc ( 'strerror', status );
		snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DIMID', ncerr );
	end
end


%
% go into define mode
status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	mexnc ( 'close', ncid );
	snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:REDEF', ncerr );
end

[varid, status] = mexnc ( 'DEF_VAR', ncid, varstruct.Name, varstruct.Nctype, num_dims, dimids );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	mexnc ( 'endef', ncid );
	mexnc ( 'close', ncid );
	snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DEF_VAR', ncerr );
end



status = mexnc ( 'enddef', ncid );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	mexnc ( 'close', ncid );
	snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:ENDDEF', ncerr );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 )
	ncerr = mexnc ( 'strerror', status );
	snc_error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:CLOSE', ncerr );
end



%
% Now just use nc_attput to put in the attributes
for j = 1:length(varstruct.Attribute)
	attname = varstruct.Attribute(j).Name;
	attval = varstruct.Attribute(j).Value;
	nc_attput ( ncfile, varstruct.Name, attname, attval );
end


return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varstruct = validate_varstruct ( varstruct )

%
% Check that required fields are there.
% Must at least have a name.
if ~isfield ( varstruct, 'Name' )
	snc_error ( 'SNCTOOLS:NC_ADDVAR:badInput', 'structure argument must have at least the ''Name'' field.' );
end

%
% Check that required fields are there.
% Default Nctype is double.
if ~isfield ( varstruct, 'Nctype' )
	varstruct.Nctype = 'double';
end

%
% Are there any unrecognized fields?
fnames = fieldnames ( varstruct );
for j = 1:length(fnames)
	fname = fnames{j};
	switch ( fname )

	case { 'Nctype', 'Name', 'Dimension', 'Attribute' }
		%
		% These are used to create the variable.  They are ok.
		
	case { 'Unlimited', 'Size', 'Rank' }
		%
		% These come from the output of nc_getvarinfo.  We don't 
		% use them, but let's not give the user a warning about
		% them either.

	otherwise
		fprintf ( 2, '%s:  unrecognized field name ''%s''.  Ignoring it...\n', mfilename, fname );
	end
end

% If the datatype is not a string.
% Change suggested by Brian Powell
if ( isa(varstruct.Nctype, 'double') && varstruct.Nctype < 7 )
	types={ 'byte' 'char' 'short' 'int' 'float' 'double'};
	varstruct.Nctype = char(types(varstruct.Nctype));
end


%
% Check that the datatype is known.
switch ( varstruct.Nctype )
case { 'NC_DOUBLE', 'double', ...
	'NC_FLOAT', 'float', ...
	'NC_INT', 'int', ...
	'NC_SHORT', 'short', ...
	'NC_BYTE', 'byte', ...
	'NC_CHAR', 'char'  }
	%
	% Do nothing
otherwise
	snc_error ( 'SNCTOOLS:NC_ADDVAR:unknownDatatype', 'unknown type ''%s''\n', mfilename, varstruct.Nctype );
end


%
% Check that required fields are there.
% Default Dimension is none.  Singleton scalar.
if ~isfield ( varstruct, 'Dimension' )
	varstruct.Dimension = [];
end

%
% Check that required fields are there.
% Default Attributes are none
if ~isfield ( varstruct, 'Attribute' )
	varstruct.Attribute = [];
end

