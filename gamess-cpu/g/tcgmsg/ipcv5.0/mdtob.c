/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv5.0/mdtob.c,v 1.1.1.5 2007-10-30 10:14:18 jmht Exp $ */

#include "sndrcv.h"

/*
  These routines use C's knowledge of the sizes of data types
  to generate a portable mechanism for FORTRAN to translate
  between bytes, integers and doubles. Note that we assume that
  FORTRAN integers are the same size as C longs.
*/

long MDTOB_(n)
     long *n;
/*
  Return the no. of bytes that n doubles occupy
*/
{
  if (*n < 0)
    Error("MDTOB_: negative argument",*n);

  return (long) (*n * sizeof(double));
}
