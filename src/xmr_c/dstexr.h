#ifndef _DSTEXR_H
#define _DSTEXR_H
//======================================================================================
//  Interface to the tridiagonal MRRR, serial version.
//
//  DSTEXR computes selected eigenpairs with indices wil:wiu for a symmetric
//  tridiagonal matrix T. The indexing of eigenpairs is with respect to an
//  ascending order of the eigenvalues, and the delivered eigenvalues will be
//  ordered ascendingly.
//
//  The computed eigenpairs will be consistent in the following sense: If the
//  routine is called twice with index sets wil1:wiu1 and wil2:wiu2, and those
//  index sets are non-overlapping, say wiu1 < wil2, then the results will obey
//   (1)  The computed eigenvalues from the first call are all <= the computed
//        eigenvalues from the second call.
//   (2)  The computed vectors from both calls will be numerically orthogonal
//        to each other.
//
//  Esp. (2) is the key to the easy parallelization approach, where each processor
//  computes one about equally sized chunk of the desired eigenpairs separately
//  from all other processors. The benefit of (1) is that the computed eigenpairs
//  won't need to be sorted across processors afterwards.
//
//  Arguments
//  =========
//  n        (input) int
//           Dimension of the matrix T, should be >= 1.
//
//  wil      (input) int
//  wiu      (input) int
//           wil:wiu is the index range of desired eigenpairs,
//           should obey 1 <= wil <= wiu <= n.
//
//  D        (input) double array, size n
//           D[1:n] are the diagonal entries of T
//
//  E        (input) double array, size n-1
//           E[1:n-1] hold the offdiagonal entries of T
//
//  W        (output) double array, dimension wiu-wil+1
//           The computed eigenvalues, in ascending order.
//           Thus, for 0 <= i < (wiu-wil+1), W[i] is an approximation to
//           eigenvalue wil+i of the matrix T.
//
//           The routine will allocate the memory and set the pointer accordingly
//           on output.
//
//  Z        (output) double array, dimension n x (wiu-wil+1)
//           The computed eigenvectors, in one chunk of column-oriented memory, that
//           is, for 0 <= i < (wiu-wil+1), the i'th computed eigenvector, which is
//           associated with W[i] and eigenvalue wil+i of the matrix, is stored in
//             Z[i*n],...,Z[(i+1)*n - 1].
//
//           The routine will allocate the memory and set the pointer accordingly
//           on output.
//
//  Return Value
//  ============
//    = -2  could not allocate enough memory
//    = -1  there was a usage error, one of the arguments had an
//          illegal value
//    =  0  means everything is ok
//    =  1  means the routine failed to converge for all eigenpairs
//          (Note: We could return status flags per eigenpair,
//           if that would be useful.)
//    >  1  some other error code
//
//--------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif 
	int dstexr_i(int n, int wil, int wiu, const double* D, const double* E, double** W, double** Z);
#ifdef __cplusplus
}
#endif 
//======================================================================================
#endif
