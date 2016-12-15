/*
 *  Main routines for the FES prediction software.
 * 
 *  File      : 
 *  Developer : CLS - CNRS/LEGOS
 *  Version   : 1.4
 *  Date      : 11 May 2005
 *  
 */

#include "fes-int.h"


/* constants */
static const int nconst = 8; 
static const int iwave[] = 
{
     0, 1, 8, 2, 4, 5, 7, 6
};



// ///////////////////////////////////////////////////////////////////////////
// Get amp and phase at a specific pos from Atlas grid for 8 constitutents 
// Data sored in table ampha that covers the hycom model boundary
// Parameters:
//   fes:	        Internal data handle, which is a pointer to a fesData
//			structure containing information about the tide
//			computation.
//   lat		Latitude in degrees (positive north) for the position
//			at which tide is computed.
//   lon		Longitude in degrees (positive north) for the position
//			at which tide is computed.
//   time		Julian day (days since 1950-01-01 00:00:00.000 UTC).
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.

int fesCore8Const(fesData* fes,
	    const double lat,
	    const double lon,
            const int ii,
            const int kk,
            const int nconst,
            const int m,
	    double    ampha[2*(ii+kk)][2*nconst])
{
    int       rc;
    int       i,k,iamp,ipha;
    double    a,theta,k1,k2,theta1,theta2;
    double    deg_conv=180/M_PI;

    
    rc = interp(fes, lat, lon);
    if ( rc != FES_SUCCESS )
	return rc;

    if ( fes->isData )
    {
         /*m_harmonic*/
        //Order constituents list for csr:
        //Q1 O1 P1 K1 N2 M2 S2 K2
        // corresponding wave code number in FES2004 is respectively
        // 0 1 8 2 4 5 7 6
        
        for ( i = 0; i < nconst; i++ )
	{

            k=iwave[i];
             // The phase and amplitudes of the constituents
            // are defined by the real and complex form of each wave
            // waves[i].cplx.re=acos(theta) &
            // waves[i].cplx.im=asin(theta)
            // we can get a and theta with selecting for a a>0 the arccos and arcsin values that matches
            k1=fes->waves[k].cplx.re;
            k2=fes->waves[k].cplx.im;
            a=sqrt(k1*k1 +k2*k2);
            theta1=acos(k1/a);
            theta2=asin(k2/a);

            if (fabs(theta1 - theta2)<EPSILON) theta=theta1;
            else if (fabs(theta1 + theta2)<EPSILON) theta=-theta1;
            else if (fabs(theta1 - M_PI + theta2)<EPSILON)  theta=theta1;  
            else if (fabs(theta1 - M_PI - theta2)<EPSILON)  theta=-theta1;
            else
            {
               printf("PB PB PB phase const!!! %f %f %f",theta,theta1,theta2);
               return -1;
            }

            if (theta<0) theta=theta + 2* M_PI;


            //update const array, in the same way than for csr
            iamp=2*i;
            ipha=2*i+1;
            ampha[m][ipha]= theta*deg_conv;       //Units mm  and dgr as it is read in mod_tides 
            ampha[m][iamp]= a*10  ;  
            //printf("%d %f %f",m,ampha[m][ipha],ampha[m][iamp]);
	}
        //printf("\n");
    }
    else
	return FES_NO_DATA;

    return FES_SUCCESS;
}


// ///////////////////////////////////////////////////////////////////////////
// Get amp and phase at a specific pos from Atlas grid for 8 constitutents 
// Data sored in table ampha that covers the hycom model boundary
// Parameters:
//   fes:	        Internal data handle, which is a pointer to a fesData
//			structure containing information about the tide
//			computation.
//   lat		Latitude in degrees (positive north) for the position
//			at which tide is computed.
//   lon		Longitude in degrees (positive north) for the position
//			at which tide is computed.
//   time		Julian day (days since 1950-01-01 00:00:00.000 UTC).
//
// Return value:
//   Returns FES_SUCCESS if the operation completed successfully or an error
//   code if a problem occurred.

int diag8Const(fesData* fes,
	    const double lat,
	    const double lon,
            const int nconst,
	    double    ampha[2*nconst])
{
    int       rc;
    int       i,k,iamp,ipha;
    double    a,theta,k1,k2,theta1,theta2;
    double    deg_conv=180/M_PI;
   
    rc = interp(fes, lat, lon);
    if ( rc != FES_SUCCESS )
	return rc;

    if ( fes->isData )
    {
         /*m_harmonic*/
        //Order constituents list for csr:
        //Q1 O1 P1 K1 N2 M2 S2 K2
        // corresponding wave code number in FES2004 is respectively
        // 0 1 8 2 4 5 7 6
       
         for ( i = 0; i < nconst; i++ )
	{
            k=iwave[i];

             // The phase and amplitudes of the constituents
            // are defined by the real and complex form of each wave
            // waves[i].cplx.re=acos(theta) &
            // waves[i].cplx.im=asin(theta)
            // we can get a and theta with selecting for a a>0 the arccos and arcsin values that matches
            k1=fes->waves[k].cplx.re;
            k2=fes->waves[k].cplx.im;
            a=sqrt(k1*k1 +k2*k2);
            theta1=acos(k1/a);
            theta2=asin(k2/a);

            if (fabs(theta1 - theta2)<EPSILON) theta=theta1;
            else if (fabs(theta1 + theta2)<EPSILON) theta=-theta1;
            else if (fabs(theta1 - M_PI + theta2)<EPSILON)  theta=theta1;  
            else if (fabs(theta1 - M_PI - theta2)<EPSILON)  theta=-theta1;
            else
            {
               printf("PB PB PB phase const!!! %f %f %f",theta,theta1,theta2);
               return -1;
            }

            if (theta<0) theta=theta + 2* M_PI;

            //update const array, in the same way than for csr
            iamp=2*i;
            ipha=2*i+1;
            ampha[ipha]= theta *deg_conv;       //Units mm in csr,dgr
            ampha[iamp]= a*10  ;  
	}
    }
    else
	return FES_NO_DATA;

    return FES_SUCCESS;
}
