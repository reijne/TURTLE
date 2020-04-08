/*
 $Id: peigs_dlasq2a.c,v 1.2 2000-10-26 15:38:25 psh Exp $
*/
#include <stdio.h>
#include <math.h>
#include "globalp.c.h"

#define SHIFT 0.0e0

void peigs_dlasq2a( Integer n, DoublePrecision *d, DoublePrecision *l, DoublePrecision *eval, DoublePrecision *work, Integer *info)
     /*
       calling routine for tridiagonal eval calculations
       in qd format
       work[10*n]
     */
{   

  Integer i, j, iii, me, nproc;
  DoublePrecision *ee, *q, *alpha, *beta, *a, *b, *zz;
  extern void dlasq2a_();
  extern void types1_();
  FILE *file;
  char filename[40];
  
  alpha = &d[0];
  beta = &l[1];
  ee = &work[0];
  q = ee + n;
  a = q + n;
  b = a + n;
  zz = b + n;
  
  me    = mxmynd_();
  nproc = mxnprc_();
  
  /*
    .... compute the q's and e's .............................
  ********/
  
  types1_(a, b, alpha, beta, q, ee, &n, &n);
  
  zz[1] = q[1];
  for ( iii = 1 ; iii <= n-1; iii++ ){
    zz[2*iii] = ee[iii];
    zz[2*iii+1] = q[iii+1];
  }
  zz[2*n] = 0.0;
  
  for ( iii = 1 ; iii <= 2*n; iii++ )
    zz[iii-1] = zz[iii];
  
  dlasq2a_( &n, &zz[0], info );
  
  if ( *info != 0 ){
    if ( me == 0 ) {
      printf(" error in dlasq2a info = %d trace %f  \n", *info, zz[2*n], zz[2*n+1] );
      sprintf( filename, "pdspevx.%d", me);
      file = fopen(filename, "w");
      for ( iii = 0; iii < n; iii++)
	fprintf(file, " %d %20.16f %20.16f \n", iii, d[iii], l[iii]);
      close(file);
      fflush(file);
      fflush(stdout);
    }
    return;
  }
  
  for(i = 0;i < n;i++)
    eval[i] = zz[n-i-1];
  
  return;
}

