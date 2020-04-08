#include <pthread.h>
#include <mpi.h>

#define TRUE  1
#define FALSE 0

pthread_t thread;     /* the thread object */
MPI_Comm thread_comm; /* communicator used to manage thread */
MPI_Status status;
int send_buffer = 0;
int recv_buffer;
int data_server_running = FALSE;
int thread_exit = 0;
int *thread_exit_ptr;
int source; /* where the terminate message will come from */
int destination; /* where the terminate message will go to */
int msg_id = 123;
int msg_len = 1;

void data_server(void)
{
   /* Provide the data server functionality by diving into and waiting in
      a MPI_recv call. Use the communicator set up in "create_thread" so
      that only the corresponding message from "destroy_thread" can 
      satisfy the MPI_recv. When the message is received it releases the
      thread leading it to terminate.
   */
   if (MPI_Send(&send_buffer,msg_len,MPI_INTEGER,source,msg_id,thread_comm))
     MPI_Abort(MPI_COMM_WORLD,999);

   pthread_exit((void*)&thread_exit);
}

void create_thread(void)
{
   /* Create a new thread to handle the MPI comms

      - create a new communicator
      - create a new thread that will just dive into and wait in MPI_recv
   */
   int colour; /* colour needed to define processor group */
   int rank = 0;

   if (data_server_running) return;
   data_server_running = TRUE;

   /* The new communicator will include only one process so set the colour
      to the rank of this process */
   if (MPI_Comm_rank(MPI_COMM_WORLD,&colour)) MPI_Abort(MPI_COMM_WORLD,905);
   if (MPI_Comm_split(MPI_COMM_WORLD,colour,rank,&thread_comm)) MPI_Abort(MPI_COMM_WORLD,901);

   source = 0;
   destination = 0;

   /* Create the thread to hang in MPI to act as a data server */
   if (pthread_create(&thread,NULL,(void *(*)(void*))data_server,NULL)) MPI_Abort(MPI_COMM_WORLD,902);
}

void destroy_thread(void)
{
   /* Destroy the thread that hangs in the MPI_recv

      - send the message that release the MPI_recv and terminates the thread
      - clean up the MPI communicator
   */

   if (!data_server_running) return;

   /* Send message to release the data server thread */
   if (MPI_Recv(&recv_buffer,msg_len,MPI_INTEGER,destination,msg_id,thread_comm,&status))
     MPI_Abort(MPI_COMM_WORLD,998);

   /* Wait for the data server thread to terminate (should not be necessary) */
   if (pthread_join(thread,(void**)&thread_exit_ptr)) MPI_Abort(MPI_COMM_WORLD,911);

   /* Tidy up the communicator we used to manage the data server thread */
   if (MPI_Comm_free(&thread_comm)) MPI_Abort(MPI_COMM_WORLD,910);

   data_server_running = FALSE;
}

/* The following instances of create_thread and destroy_thread are
   simply aliases to address linking issues with Fortran codes.
   In this case that is the simplest way to handle this as there are
   no arguments anyway.
*/

void create_thread_(void)
{
    create_thread();
}
void destroy_thread_(void)
{
    destroy_thread();
}
 
void create_thread__(void)
{
    create_thread();
}
void destroy_thread__(void)
{
    destroy_thread();
}
