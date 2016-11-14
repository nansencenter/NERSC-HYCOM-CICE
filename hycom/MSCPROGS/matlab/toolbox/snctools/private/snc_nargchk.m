function snc_nargchk(low,high,N)
% SNC_NARGCHK:  wrapper for NARGCHK, which changed functionality at R???

v = version;

switch ( version('-release') )
	case { '12', '13' }
		error ( nargchk(low,high,N ) );
	otherwise
		error ( nargchk(low,high,N,'struct') );
end


