/*
When building the GA/MPI/ScaLAPACK code on BlueGene the following arises:
1. GAMESS-UK has to be build with -qextname otherwise the Global Arrays will
   not link.
2. The BlueGene MPI libraries do not provide interfaces with trailing
   underscores.
The way out of this problem is provided by the functions in this file
which do nothing but map MPI routines with trailing underscores to MPI
routines without trailing underscores. These routines are written in C to 
avoid Fortran compile flags messing with what we are trying to achieve here.
All very boring but effective.
*/

double mpi_wtime(void);

void mpi_init_(int *ierr)
{
   mpi_init(ierr);
}

void mpi_initialized_(int* flag, int *ierr)
{
   mpi_initialized(flag,ierr);
}

void mpi_finalize_(int *ierr)
{
   mpi_finalize(ierr);
}

void mpi_abort_(int *comm, int *type, int *ierr)
{
   mpi_abort(comm,type,ierr);
}

void mpi_attr_get_(int *comm, int *keyval, int *attrval, int *flag, int *ierr)
{
   mpi_attr_get(comm,keyval,attrval,flag,ierr);
}

void mpi_barrier_(int* comm, int *ierr)
{
   mpi_barrier(comm,ierr);
}

void mpi_bcast_(void *buf, int *count, int *datatype, int *root, int *comm,
                int *ierr)
{
   mpi_bcast(buf,count,datatype,root,comm,ierr);
}

void mpi_comm_create_(int *comm, int *group, int *newgroup, int *ierr)
{
   mpi_comm_create(comm,group,newgroup,ierr);
}

void mpi_comm_dup_(int *comm, int *newcomm, int *ierr)
{
   mpi_comm_dup(comm,newcomm,ierr);
}

void mpi_comm_free_(int *comm, int *ierr)
{
   mpi_comm_free(comm,ierr);
}

void mpi_comm_group_(int *comm, int *group, int *ierr)
{
   mpi_comm_group(comm,group,ierr);
}

void mpi_comm_rank_(int *comm, int *rank, int *ierr)
{
   mpi_comm_rank(comm,rank,ierr);
}

void mpi_comm_size_(int *comm, int *size, int *ierr)
{
   mpi_comm_size(comm,size,ierr);
}

void mpi_comm_split_(int *comm, int *color, int *key, int *newcomm, int *ierr)
{
   mpi_comm_split(comm,color,key,newcomm,ierr);
}

void mpi_error_class_(int *errorcode, int *errorclass, int *ierr)
{
   mpi_error_class(errorcode,errorclass,ierr);
}

void mpi_group_free_(int *group, int *ierr)
{
   mpi_group_free(group,ierr);
}

void mpi_group_incl_(int *group, int *n, int *ranks, int *newgrp, int *ierr)
{
   mpi_group_incl(group,n,ranks,newgrp,ierr);
}

void mpi_group_translate_ranks_(int *grp1, int *n, int *ranks1, int *grp2, 
                                int *ranks2, int *ierr)
{
   mpi_group_translate_ranks(grp1,n,ranks1,grp2,ranks2,ierr);
}

void mpi_op_create_(void *function, int *commute, int *op, int *ierr)
{
   mpi_op_create(function,commute,op,ierr);
}

void mpi_op_free_(int *op, int *ierr)
{
   mpi_op_free(op,ierr);
}

void mpi_pack_(void *inbuf, int *incount, int *datatype, void *outbuf,
               int *outcount, int *position, int *comm, int *ierr)
{
   mpi_pack(inbuf,incount,datatype,outbuf,outcount,position,comm,ierr);
}

void mpi_pack_size_(int *incount, int *datatype, int *comm, int *size,
                    int *ierr)
{
   mpi_pack_size(incount,datatype,comm,size,ierr);
}

void mpi_send_(void *buf, int *count, int *datatype, int *dest, int *tag, 
               int *comm, int *ierr)
{
   mpi_send(buf,count,datatype,dest,tag,comm,ierr);
}

void mpi_rsend_(void *buf, int *count, int *datatype, int *dest, int *tag,
                int *comm, int *ierr)
{
   mpi_rsend(buf,count,datatype,dest,tag,comm,ierr);
}

void mpi_sendrecv_(void *sendbuf, int *sendcount, int *sendtype, int *dest,
                   int  *sendtag,
                   void *recvbuf, int *recvcount, int *recvtype, int *source,
                   int  *recvtag,
                   int  *comm, int *ierr)
{
   mpi_sendrecv(sendbuf,sendcount,sendtype,dest,sendtag,
             recvbuf,recvcount,recvtype,source,recvtag,
             comm,ierr);
}

void mpi_isend_(void *buf, int *count, int *datatype, int *dest, int *tag,
                int *comm, int *request, int *ierr)
{
   mpi_isend(buf,count,datatype,dest,tag,comm,request,ierr);
}

void mpi_recv_(void *buf, int *count, int *datatype, int *source, int *tag,
              int *comm, int *status, int *ierr)
{
   mpi_recv(buf,count,datatype,source,tag,comm,status,ierr);
}

void mpi_irecv_(void *buf, int *count, int *datatype, int *dest, int *tag,
                int *comm, int *request, int *ierr)
{
   mpi_irecv(buf,count,datatype,dest,tag,comm,request,ierr);
}

void mpi_allreduce_(void *sendbuf, void *recvbuf, int *count, int *datatype,
                    int  *op, int *comm, int *ierr)
{
   mpi_allreduce(sendbuf,recvbuf,count,datatype,op,comm,ierr);
}

void mpi_reduce_(void *sendbuf, void *recvbuf, int *count, int *datatype,
                 int  *op, int *root, int *comm, int *ierr)
{
   mpi_reduce(sendbuf,recvbuf,count,datatype,op,root,comm,ierr);
}

void mpi_testall_(int *count, int *array_of_requests, int *flag,
                  int *array_of_statuses, int *ierr)
{
   mpi_testall(count,array_of_requests,flag,array_of_statuses,ierr);
}

void mpi_type_commit_(int *type, int *ierr)
{
   mpi_type_commit(type,ierr);
}

void mpi_type_free_(int *datatype, int *ierr)
{
   mpi_type_free(datatype,ierr);
}

void mpi_type_indexed_(int *count, int *blocklens, int *indices, int *old_type,
                       int *newtype, int *ierr)
{
   mpi_type_indexed(count,blocklens,indices,old_type,newtype,ierr);
}

void mpi_type_struct_(int *count, int *blocklens, int *indices, int *old_type,
                      int *newtype, int *ierr)
{
   mpi_type_struct(count,blocklens,indices,old_type,newtype,ierr);
}

void mpi_type_vector_(int *count, int *blocklength, int *stride, int *old_type,
                      int *newtype, int *ierr)
{
   mpi_type_vector(count,blocklength,stride,old_type,newtype,ierr);
}

void mpi_wait_(int *request, int *status, int *ierr)
{
   mpi_wait(request,status,ierr);
}

void mpi_waitall_(int *count, int *array_of_requests,
                  int *array_of_statuses, int *ierr)
{
   mpi_waitall(count,array_of_requests,array_of_statuses,ierr);
}


double mpi_wtime_(void) 
{
   return mpi_wtime();
}
