#
# --- awk script that converts HYCOM archive names from new to old format.
#

	{
	y = substr( $0,  7, 4 )
	d = substr( $0, 12, 3 )
	m = (y-1)*366 + d - 16
	printf("archk.%6.6d\n",m)
}
