/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv5.0/drand48.c,v 1.1.1.5 2007-10-30 10:14:17 jmht Exp $ */

#include "srftoc.h"

extern long random();
extern int srandom();

double DRAND48_()
{
  return ( (double) random() ) * 4.6566128752458e-10;
}

void SRAND48_(seed)
  unsigned *seed;
{
  (void) srandom(*seed);
}
