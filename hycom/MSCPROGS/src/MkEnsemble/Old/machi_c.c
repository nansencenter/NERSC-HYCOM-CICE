#if defined(SGI)
/* 
 --- Fortran-callable routine ZUNDER that sets the bit to specify
 --- that underflows are flushed to zero in hardware on SGI R10000.
 --- See man handle_sigfpes
 --- Alan J. Wallcraft,  NRL,  October 1997.
*/
#include <sys/fpu.h>
void zunder_()
{
	union fpc_csr   n;
	n.fc_word = get_fpc_csr();
	n.fc_struct.flush = 1;
	set_fpc_csr(n.fc_word);
}
#endif /* SGI */

#if defined(AIX)
/* 
 --- Fortran-callable function WTIME that returns the wall time in seconds.
 --- Probably not thread-safe, only for Power-PC systems.
 --- Alan J. Wallcraft,  NRL,  May 2001.
 --- Based on notes by Bob Walkup (10x faster than MPI_WTIME).
*/
#include <sys/time.h>
#include <sys/systemcfg.h>
double wtime(void)
{
	struct timebasestruct TB;
	static int    first_call;
	static double tb_factor;
	double tb_top,tb_bot;
	if (first_call == 0) {
		first_call = 1;
		tb_top     = (double) _system_configuration.Xint;
		tb_bot     = (double) _system_configuration.Xfrac;
		tb_factor  = tb_top/tb_bot;
        }
	read_real_time(&TB, TIMEBASE_SZ);
	return ( tb_factor * ( 4.294967296*((double) TB.tb_high) + 1.0e-9*((double) TB.tb_low) ) );
}
#endif /* AIX */

void machine_c()
{
}
