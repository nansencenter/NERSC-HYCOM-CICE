function validate_index_vectors(start,count,stride,nvdims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% VALIDATE_INDEX_VECTORS
%    We need to check the lengths of the index vectors before calling the
%    netCDF library.  A bad length can confuse the mex-file, and it's really 
%    not the mex-file's responsibility to do this.  We can't do this when 
%    parsing the input arguments since we need the number of dimensions 
%    corresponding to the input variable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(start) && (length(start) ~= nvdims)
    mexnc('close',ncid);
    fmt = 'length of the start index vector (%d) does not equal the number of dimensions (%d) for %s';
    msg = sprintf ( fmt, length(start), nvdims, varname );
    snc_error ( 'SNCTOOLS:validate_index_vectors:badStartIndex', msg );
end
if ~isempty(count) && (length(count) ~= nvdims)
    mexnc('close',ncid);
    fmt = 'length of the count index vector (%d) does not equal the number of dimensions (%d) for %s';
    msg = sprintf ( fmt, length(count), nvdims, varname );
    error ( 'SNCTOOLS:validate_index_vectors:badCountIndex', msg );
end
if ~isempty(stride) && (length(stride) ~= nvdims)
    mexnc('close',ncid);
    fmt = 'length of the stride index vector (%d) does not equal the number of dimensions (%d) for %s';
    msg = sprintf ( fmt, length(count), nvdims, varname );
    error ( 'SNCTOOLS:validate_index_vectors:badStrideIndex', msg );
end

%
% Don't bother if any of COUNT == 0
if prod(count) == 0
    snc_error ( 'SNCTOOLS:validate_index_vectors:badCountIndex', 'retrieval request was for 0 elements.' );  
end


return






