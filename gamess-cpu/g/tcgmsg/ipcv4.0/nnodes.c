/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/nnodes.c,v 1.1.1.5 2007-10-30 10:14:14 jmht Exp $ */

#include "sndrcv.h"
#include "sndrcvP.h"

long NNODES_()
/*
  return total no. of processes
*/
{
  return SR_n_proc;
}

