#include "dstexr.h"

#include "c2fort.h"

#include <stdlib.h>

//#include <new>
//======================================================================================
// Declare interface to the Fortran routine
#ifdef __cplusplus
extern "C" {
#endif
	void
	EXTFORTNAME(dstexr)(
        int* n, double* D, double* E, int* wil, int* wiu,
        double* W, double* Z, int* ldz, int* ISuppZ,
        double* RWork, int* lrwork, int* IWork, int* liwork,
        int* info
    );
#ifdef __cplusplus
}
#endif
//======================================================================================
int
dstexr_i(
        int n, int wil, int wiu,
        const double* D, const double* E,
        double** W, double** Z
    )
{
    //-- Check Arguments
    if( n <= 1 || !(1 <= wil && wil <= wiu && wiu <= n) ||
        D == 0 || E == 0 )
    {
        return -1;
    }

    //-- Workspace Query
    int liwork = -1, lrwork = -1;

    int idum = -1; double rdum = 0.0;
    EXTFORTNAME(dstexr)(
        &n, &rdum, &rdum, &idum, &idum, &rdum, &rdum, &idum, &idum,
        &rdum, &lrwork, &idum, &liwork, &idum );

    //-- Allocate Memory
    const int wlen = wiu - wil + 1;
    int res = 0;
    double* RWork = 0;
    int* IWork = 0;
    int* ISuppZ = 0;
    *W = 0;
    *Z = 0;

    *Z = (double*)malloc(n * wlen * sizeof(double));
    if(*Z == 0){free(*Z); return -2;}

    *W = (double*)malloc(wlen * sizeof(double));
    if(*W == 0){ return -2;}

    RWork = (double*)malloc(lrwork * sizeof(double));
    IWork = (int*)malloc(liwork * sizeof(int));

    ISuppZ = (int*)malloc(2 * wlen * sizeof(int));
    if(*Z == 0 || *W == 0 || RWork == 0 || IWork == 0 || ISuppZ == 0){
		free(*W);
		free(RWork);
		free(IWork);
		free(ISuppZ); return -2;
	}

    //-- The Call
    int info;

    double* pD = (double*)D;
    double* pE = (double*)E;

    EXTFORTNAME(dstexr)(
            &n, pD, pE, &wil, &wiu, *W, *Z, &n, ISuppZ,
            RWork, &lrwork, IWork, &liwork, &info );

    return info;
}
//======================================================================================
