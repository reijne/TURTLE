_IF(parallel)
#include <mpi.h>

_IF(tcgmsg-mpi)
#include "tcgmsgP.h"
_ENDIF

/* very simple mpi_sendrecv to link in with GA jvl,2005 */
void sendrecv_(stype, sbuf, lsbuf, snode, rtype, rbuf, lrbuf, rnode)
     int  *stype;
_IF(tcgmsg-mpi)
     Void *sbuf;
_ELSEIF(_AND(ga,mpi))
     void *sbuf;
_ENDIF
     int  *lsbuf;
     int  *snode;
     int  *rtype;
_IF(tcgmsg-mpi)
     Void *rbuf;
_ELSEIF(_AND(ga,mpi))
     void *rbuf;
_ENDIF
     int  *lrbuf;
     int  *rnode;
{
  int ierr;
  MPI_Status status;
  /* For the time being, we assume that we are being compiled with 
     i4 and so that an int is suitable for holding the communicator */
  MPI_Comm Comm;
_IF(tcgmsg-mpi)
  Comm = TCGMSG_Comm;
_ELSEIF(_AND(ga,mpi))
  extern void getcomms_(int *cworld, int *cgamess, int * cworkers);
  int commworld, commgamess, commworkers;
  getcomms_(&commworld, &commgamess, &commworkers);
  Comm = commgamess;
_ENDIF

    ierr = MPI_Sendrecv(sbuf, (int)*lsbuf, MPI_CHAR, (int)*snode, (int)*stype,
		        rbuf, (int)*lrbuf, MPI_CHAR, (int)*rnode, (int)*rtype, 
			Comm, &status);

_IF(tcgmsg-mpi)
    tcgmsg_test_statusM("SENDRECV_:", ierr);
_ENDIF
}
_ELSE
void func_with_no_name()
{
            return;
}
_ENDIF
