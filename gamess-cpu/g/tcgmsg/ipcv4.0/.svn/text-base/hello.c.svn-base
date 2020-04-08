/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/hello.c,v 1.1.1.5 2007-10-30 10:14:14 jmht Exp $ */

#include "sndrcv.h"

int main(argc, argv)
     int argc;
     char **argv;
/*
  Traditional first parallel program
*/
{
  PBEGIN_(argc, argv);

  (void) printf("Hello from node %ld\n",NODEID_());

  PEND_();

  return 0;
}
