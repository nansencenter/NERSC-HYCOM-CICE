function snc_error ( error_id, error_msg )

v = version;

switch ( version('-release') )
	case '12'
		error ( error_msg );
	otherwise
		error ( error_id, error_msg );
end
