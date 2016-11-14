function theBuffer = nc_getbuffer ( ncfile, varargin )
% NC_GETBUFFER:  read the unlimited variables of a netcdf file into a structure
%
% USAGE:  theBuffer = nc_getbuffer ( ncfile );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, varlist );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, start, count );
% USAGE:  theBuffer = nc_getbuffer ( ncfile, varlist, start, count );
%
% PARAMETERS:
% INPUT:
%     ncfile:  
%        Input netcdf file name.
%     varlist:
%        cell array of named variables.  Only data for these variables 
%        will be retrieved from the file
%     start, count:
%        starting index and number of records to retrieve.  This is 
%        optional.  If not provided, all of the record variables will
%        be retrieved.
%        
%        If start is negative, then the last few records (total of
%        "count") are retrieved.
%
%        If count is negative, then everything beginning at "start" 
%        and going all the way to the end is retrieved.
% 
% OUTPUT:
%     theBuffer:
%        Structure with fields corresponding to each netcdf record variable.  
%        Each such field contains the data for that variable.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_getbuffer.m 2315 2007-09-03 16:07:33Z johnevans007 $
% $LastChangedDate: 2007-09-03 12:07:33 -0400 (Mon, 03 Sep 2007) $
% $LastChangedRevision: 2315 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
% assume failure until success is known
theBuffer = [];

snc_nargchk(1,4,nargin);
snc_nargoutchk(1,1,nargout);


%
% check that the first argument is a char
if ~ischar ( ncfile )
   	snc_error (  'SNCTOOLS:NC_GETBUFFER:badInput', 'filename argument must be character.' );
end


[varlist,start,count] = parse_inputs ( varargin{:} );

metadata = nc_info ( ncfile );

num_datasets = length(metadata.Dataset);

skip_this_variable = construct_skip_list(varlist,metadata);

%
% Find the unlimited dimension and it's length
record_length = -1;
num_dims = length(metadata.Dimension);
for j = 1:num_dims
	if metadata.Dimension(j).Unlimited
		record_length = metadata.Dimension(j).Length;
	end
end
if record_length < 0
   	snc_error (  'SNCTOOLS:NC_GETBUFFER:noUnlimitedDimension', ...
	         'An unlimited dimension is required.');
end

%
% figure out what the start and count really are.
if ~isempty(start) && ~isempty(count)
	if start < 0
		start = record_length - count;
	end
	if count < 0
		count = record_length - start;
	end
	if (start < 0) && (count < 0)
   		snc_error (  'SNCTOOLS:NC_GETBUFFER:badIndexing', ...
	             'both start and count cannot be less than zero.');
	end
end




for j = 1:num_datasets
	
	%
	% Did we restrict retrieval to a few variables?
	if ~isempty(varlist) && skip_this_variable(j)
		continue
	end

	%
	% If it is not an unlimited variable, we don't want it.
	if ~metadata.Dataset(j).Unlimited
		continue
	end

	if ~isempty(start) && ~isempty(count) 
		varstart = zeros(size(metadata.Dataset(j).Size));
		varstart(1) = start;
		varcount = metadata.Dataset(j).Size;
		varcount(1) = count;

		vardata = nc_varget ( ncfile, metadata.Dataset(j).Name, varstart, varcount );

	else
		vardata = nc_varget ( ncfile, metadata.Dataset(j).Name );
	end


	theBuffer.(metadata.Dataset(j).Name) = vardata;

end


return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varlist, start, count] = parse_inputs ( varargin )

varlist = {};
start = [];
count = [];

%
% figure out what the inputs actually were
switch nargin
case 1
	if iscell(varargin{1})
		varlist = varargin{1};
	else
		snc_error ( 'SNCTOOLS:NC_GETBUFFER:badInput', '2nd of two input arguments must be a cell array.' );
	end
case 2
	if isnumeric(varargin{1}) && isnumeric(varargin{2})
		start = varargin{1};
		count = varargin{2};
	else
		snc_error ( 'SNCTOOLS:NC_GETBUFFER:badInput', '2nd and 3rd of three input arguments must be numeric.' );
	end
case 3
	if iscell(varargin{1})
		varlist = varargin{1};
	else
		snc_error ( 'SNCTOOLS:NC_GETBUFFER:badInput', '2nd of four input arguments must be a cell array.' );
	end
	if isnumeric(varargin{2}) && isnumeric(varargin{3})
		start = varargin{2};
		count = varargin{3};
	else
		snc_error ( 'SNCTOOLS:NC_GETBUFFER:badInput', '3rd and 4th of four input arguments must be numeric.' );
	end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function skip_it = construct_skip_list(varlist,metadata)

num_datasets = length(metadata.Dataset);
if ~isempty(varlist) 

	skip_it = ones(num_datasets,1);

	%
	% Go thru and quickly set up a flag for each Dataset
	for j = 1:num_datasets
		for k = 1:length(varlist)
			if strcmp(varlist{k}, metadata.Dataset(j).Name)
				skip_it(j) = 0;
			end
		end
	end

else
	skip_it = zeros(num_datasets,1);
end


retrievable_datasets = find(1 - skip_it);
if ~any(retrievable_datasets)
	error ( 'No datasets found.\n' );
end


