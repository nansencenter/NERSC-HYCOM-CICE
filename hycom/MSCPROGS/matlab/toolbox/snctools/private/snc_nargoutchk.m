function out = snc_nargoutchk(low,high,N)
% SNC_NARGOUTCHK:  wrapper for NARGOUTCHK, which changed functionality at R???

v = version;

switch ( version('-release') )
	case { '12', '13' }
		error ( nargoutchk(low,high,N ) );
	otherwise
		error ( nargoutchk(low,high,N,'struct') );
end



