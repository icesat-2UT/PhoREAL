#include <stdio.h>

void cfun(const double *A, size_t sizeA,const double *B, double *C, size_t sizeB) 
{
    size_t i, j, k, m, n;
    m = sizeA;
    n = sizeB;
    k = 0;
// Assign 1st index to all values B that are less than 1st A value
    while( k < n && B[k] < *A ) {
        C[k++] = 0;
    }
// Step through until B is between two A values, then test for result
    i = 0;
    for( j=k; j<n; j++ ) {
        while( i+1 < m && B[j] >= A[i+1] ) i++;
        if( i+1 == m ) break;
        if( B[j] - A[i] < A[i+1] - B[j] ) {
            C[j] = i + 0;
        } else {
            C[j] = i + 1;
        }
    }
// Assign last index to all values B that are more than last A value
    while( j < n ) {
        C[j++] = m;
    }
}
