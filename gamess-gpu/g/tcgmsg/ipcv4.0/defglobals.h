/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/defglobals.h,v 1.1.1.6 2007-10-30 10:14:13 jmht Exp $ */

#ifndef SNDRCVP
#include "sndrcvP.h"
#endif

/* Actual definition of these globals ... need this once in
   any executable ... included by cluster.c */

/*********************************************************
  Global information and structures ... all begin with SR_
  ********************************************************/

long SR_n_clus;                   /* No. of clusters */
long SR_n_proc;                   /* No. of processes excluding dummy
				     master process */

int  SR_socks[MAX_PROCESS];
int  SR_socks_proc[MAX_PROCESS];
int  SR_nsock;
long SR_using_shmem;

long SR_clus_id;                  /* Logical id of current cluster */
long SR_proc_id;                  /* Logical id of current process */

long SR_debug;                    /* flag for debug output */

long SR_parallel;		/* True if job started with parallel */

long SR_exit_on_error;            /* flag to exit on error */
long SR_error;                    /* flag indicating error has been called
                                     with SR_exit_on_error == FALSE */

long SR_numchild;                   /* no. of forked processes */
long SR_pids[MAX_SLAVE];          /* pids of forked processes */


/* This is used to store info from the PROCGRP file about each
   cluster of processes */

struct cluster_info_struct SR_clus_info[MAX_CLUSTER];

struct process_info_struct SR_proc_info[MAX_PROCESS];
