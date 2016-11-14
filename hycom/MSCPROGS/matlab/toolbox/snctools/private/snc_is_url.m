function b = snc_is_url ( urlCandidate )
% SNC_IS_URL:  determines whether or not something is a URL.
%
% USAGE:  b = snc_is_url ( urlCandidate );
%
% PARAMETERS:
% Input:
%     urlCandidate:
%         We want to know if this is a URL or not.
% Output:
%     b:
%         0 if no, 1 if yes.


% 
% b == 0:  not a url
b = isfinite ( regexp ( urlCandidate, 'http://.*/.*' ) );
return
