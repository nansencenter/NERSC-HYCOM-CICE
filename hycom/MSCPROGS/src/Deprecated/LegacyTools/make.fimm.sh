OPTS="-convert big_endian -real_size 64 -double_size 64" 


ifort -o read_average $OPTS mod_dimensions.F90 mod_average.F90 p_read_average.F90 

#ifort -o datetojul m_datetojulian.F90 p_datetojul.F90 
#rm m_datetojulian.mod

