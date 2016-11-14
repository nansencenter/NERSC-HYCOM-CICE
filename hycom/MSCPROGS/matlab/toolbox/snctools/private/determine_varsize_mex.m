function the_var_size = determine_varsize_mex ( ncid, dimids, nvdims )
% DETERMINE_VARSIZE_MEX: Need to figure out just how big the variable is.
%
% VAR_SIZE = DETERMINE_VARSIZE_MEX(NCID,DIMIDS,NVDIMS);

%
% If not a singleton, we need to figure out how big the variable is.
if nvdims == 0
    the_var_size = 1;
else
    the_var_size = zeros(1,nvdims);
    for j=1:nvdims,
        dimid = dimids(j);
        [dim_size,status]=mexnc('inq_dimlen', ncid, dimid);
        if ( status ~= 0 )
            ncerr = mexnc ( 'strerror', status );
            mexnc('close',ncid);
            error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_DIM_LEN', ncerr );
        end
        the_var_size(j)=dim_size;
    end
end

return





