/* $$  */
#include "tcgmsgP.h"

static const long false = 0;
static const long true  = 1;

extern void Busy(int);
extern void flush_send_q(void);

#ifdef SHMEM

/* All data movement to/from shared memory is done using the
   COPY_TO/FROM_SHMEM macros */

#ifdef KSR
       extern copyto(const void *,  void *, long);
       extern void copyfrom(const void *,  void *, long);
#      define COPY_TO_REMOTE(src, dest, n, destnode) copyto(src, dest, n)
#      define COPY_FROM_REMOTE(src, dest, n, fromnode) copyfrom(src, dest, n)
#      define COPY_TO_LOCAL(src, dest, n) copyto(src, dest, n)
#      define COPY_FROM_LOCAL(src, dest, n) copyfrom(src, dest, n)
#elif defined(CRAY_T3E)
#      include <mpp/shmem.h>
       /* T3E has adaptive routing, this is safer with shmem_quiet */
#      define COPY_TO_REMOTE(src, dest, n, node){ \
              shmem_put((long*)(dest),(long*)(src),(int) ((n)>>3),(node));\
              shmem_quiet(); \
              }

#      define COPY_FROM_REMOTE(src, dest, n, node) \
              shmem_get((long*)(dest),(long*)(src),(int) ((n)>>3),(node))

#      define COPY_TO_LOCAL(src, dest, n)\
              (void)memcpy(dest, src, (size_t) n)

#      define COPY_FROM_LOCAL(src, dest, n)\
              (void)memcpy(dest, src, (size_t) n)

#elif defined(CRAY_T3D)
#      include <mpp/shmem.h>
#      ifdef FLUSHCACHE
#            define FLUSH_CACHE           shmem_udcflush()
#            define FLUSH_CACHE_LINE(x)   shmem_udcflush_line((long*)(x))
#      else
#            define FLUSH_CACHE 
#            define FLUSH_CACHE_LINE(x)  
#      endif

#      define COPY_TO_REMOTE(src, dest, n, node){ \
              shmem_put((long*)(dest),(long*)(src),(int) ((n)>>3),(node));\
              }

#      define COPY_FROM_REMOTE(src, dest, n, node) \
              shmem_get((long*)(dest),(long*)(src),(int) ((n)>>3),(node))

#      define COPY_TO_LOCAL(src, dest, n)\
              (void)memcpy(dest, src, (size_t) n)
/*              (void) copyto(src, dest, (long) n)*/

#      define COPY_FROM_LOCAL(src, dest, n)\
              (void)memcpy(dest, src, (size_t) n)
/*              (void) copyfrom(src, dest, (long) n)*/

#else
#define COPY_TO_LOCAL(src, dest, n, destnode) (void) memcpy(dest, src, n)
#define COPY_FROM_LOCAL(src, dest, n) (void)memcpy(dest, src, n)
#define COPY_FROM_REMOTE(src,dest,n,p) (void)memcpy(dest, src, n) 
#define COPY_TO_REMOTE(src,dest,n,p) (void)memcpy(dest, src, n)
#endif
#ifndef FLUSH_CACHE 
#	define FLUSH_CACHE          
#endif
#ifndef FLUSH_CACHE_LINE
#     define FLUSH_CACHE_LINE(x) 
#endif
#endif

static long remote_flag(long *p, long node)
/*
  Return the value of a volatile variable in shared memory
  that is REMOTE to this processor
*/
{
  long tmp;

/*  FLUSH_CACHE;*/                       /* no need to flush for one word only*/
  COPY_FROM_REMOTE(p, &tmp, sizeof(tmp), node);
  return tmp;
}

static long local_flag(long *p)
/*
  Return the value of a volatile variable in shared memory
  that is LOCAL to this processor
*/
{
  FLUSH_CACHE_LINE(p);  
  return(*p);
}

static void local_await(long *p, long value)
/*
  Wait for (*p == value)
*/
{
  long pval;
  long nspin = 0;
#ifdef NOSPIN
  long spinlim = 100;
# ifdef CRAY
  long waittim = 10000;
# endif
#else
  long spinlim = 100000000;
# ifdef CRAY
  long waittim = 100000;
# endif
#endif  

  while ((pval = local_flag(p)) != value) {

    if (pval && (pval != value)) {
      fprintf(stdout,"%2ld: invalid value=%ld, local_flag=%p %ld\n", 
              TCGMSG_nodeid, value, p, pval);
      fflush(stdout);
      exit(1);
    }
    nspin++;
    if((nspin&7)==0)flush_send_q();
    if (nspin < spinlim)
      Busy(100);
    else 
#ifdef CRAY
      USleep(waittim);
#else
      usleep(1);
#endif
  }
}

#define ABS(a) (((a) >= 0) ? (a) : (-(a)))

long async_send(SendQEntry *entry)
/*
  Entry points to info about a message ... determine which
  transport mechanism is appropriate and send as much as
  possible without blocking.

  Right now just shared memory ... when sockets are working this
  routine will become async_shmem_send.

  Shared-memory protocol aims for low latency. Each process has
  one buffer for every other process.  Thus, to send a message U
  merely have to determine if the receivers buffer for you is empty
  and copy directly into the receivers buffer.

  Return 0 if more data is to be sent, 1 if the send is complete.
*/
{
  long node = entry->node;
  ShmemBuf *sendbuf= TCGMSG_proc_info[node].sendbuf;
  long nleft, ncopy;
  long pval;
  long info[4];
  
#ifdef DEBUG
  (void) fprintf(stdout,"%2ld: sending to %ld buf=%lx len=%ld\n",
                 TCGMSG_nodeid, node, entry->buf, entry->lenbuf); 
  (void) fprintf(stdout,"%2ld: sendbuf=%lx\n", TCGMSG_nodeid, sendbuf);
  (void) fflush(stdout);
#endif

  if ((pval = remote_flag(&sendbuf->info[3], node))) {
#ifdef DEBUG
  {
    long info[4];
    FLUSH_CACHE;
    COPY_FROM_REMOTE(sendbuf->info, info, sizeof(info), node);
    fprintf(stdout,"%2ld: snd info after full = %ld %ld %ld\n", TCGMSG_nodeid,
            info[0], info[1], info[2]);
    fflush(stdout);
  }
    sleep(1);
#endif

    return 0;
  }

  info[0] = entry->type; info[1] = entry->lenbuf; info[2] = entry->tag;

  /* Copy over the first buffer load of the message */

  nleft = entry->lenbuf - entry->written;
  ncopy = (long) ((nleft <= SHMEM_BUF_SIZE) ? nleft : SHMEM_BUF_SIZE);
  
  if (ncopy&7) {
#ifdef DEBUG
    printf("%2ld: rounding buffer up %ld->%ld\n", TCGMSG_nodeid,
           ncopy, ncopy + 8 - (ncopy&7));
    fflush(stdout);
#endif
    ncopy = ncopy + 8 - (ncopy&7); 
  }

  if (ncopy) {
    COPY_TO_REMOTE(entry->buf+entry->written, sendbuf->buf, ncopy, node);
  }


  /* NOTE that SHMEM_BUF_SIZE is a multiple of 8 by construction so that
     this ncopy is only rounded up on the last write */

  ncopy = (long) ((nleft <= SHMEM_BUF_SIZE) ? nleft : SHMEM_BUF_SIZE);
  entry->written += ncopy;
  entry->buffer_number++;
  
  /* Copy over the header information include buffer full flag */

  info[3] = entry->buffer_number;
#ifdef CRAY_T3E
  COPY_TO_REMOTE(info, sendbuf->info, 3*sizeof(long), node);
  shmem_long_p(&sendbuf->info[3],info[3], node);
  shmem_quiet();
#else
  COPY_TO_REMOTE(info, sendbuf->info, sizeof(info), node);
#endif

  return (long) (entry->written == entry->lenbuf);
}

void msg_rcv(long type, char *buf, long lenbuf, long *lenmes, long node)
/*
  Receive a message of given type from the specified node, returning
  the message and length of the message.
  
  Right now just shared memory ... when sockets are working this
  routine will become msg_shmem_rcv

  Shared-memory protocol aims for low latency. Each process has
  one buffer for every other process.  Thus, to send a message U
  merely have to determine if the receivers buffer for you is empty
  and copy directly into the receivers buffer.

  Return 0 if more data is to be sent, 1 if the send is complete.
*/
{
  long me = TCGMSG_nodeid;
  ShmemBuf *recvbuf;		/* Points to receving buffer */
  long nleft;
  long msg_type, msg_tag, msg_len;
  long buffer_number = 1;
  long expected_tag = TCGMSG_proc_info[node].tag_rcv++;

  if (node<0 || node>=TCGMSG_nnodes)
    Error("msg_rcv: node is out of range", node);

  recvbuf = TCGMSG_proc_info[node].recvbuf;  

  /* Wait for first part message to be written */

#ifdef DEBUG
  (void) fprintf(stdout,"%2ld: receiving from %ld buf=%lx len=%ld\n",
                 me, node, buf, lenbuf); 

  (void) fprintf(stdout,"%2ld: recvbuf=%lx\n", me, recvbuf);
  (void) fflush(stdout);
#endif

  local_await(&recvbuf->info[3], buffer_number);

  /* Copy over the header information */

/*  FLUSH_CACHE;*/
  msg_type = recvbuf->info[0]; 
  msg_len  = recvbuf->info[1];
  msg_tag  = recvbuf->info[2];

  /* Check type and size information */

  if (msg_tag != expected_tag) {
    (void) fprintf(stdout,
		   "rcv: me=%ld from=%ld type=%ld, tag=%ld, expectedtag=%ld\n",
		   me, node, type, msg_tag, expected_tag);
    fflush(stdout);
    Error("msg_rcv: tag mismatch ... transport layer failed????", 0L);
  }

  if (msg_type != type) {
    (void) fprintf(stdout,
		   "rcv: me=%ld from=%ld type=(%ld != %ld) tag=%ld len=%ld\n",
                   me, node, type, msg_type, msg_tag, msg_len);
    fflush(stdout);
    Error("msg_rcv: type mismatch ... strong typing enforced\n", 0L);
  }

  if (msg_len > lenbuf) {
    (void) fprintf(stderr,
		   "rcv: me=%ld from=%ld type=%ld tag=%ld len=(%ld > %ld)\n",
		   me, node, type, msg_tag, msg_len, lenbuf);
    Error("msg_rcv: message too long for buffer\n", 0L);
  }
    
  nleft = *lenmes = msg_len;
  if (nleft == 0) {
   recvbuf->info[3] = false;
  }

  while (nleft) {
    long ncopy = (long) ((nleft <= SHMEM_BUF_SIZE) ? nleft : SHMEM_BUF_SIZE);
    { 
    long line;
    if(ncopy < 321) 
        for(line = 0; line < ncopy; line+=32) 
            FLUSH_CACHE_LINE(recvbuf->buf+line);
    else 
      FLUSH_CACHE;
    }

/*    if (buffer_number > 1) FLUSH_CACHE;*/
    COPY_FROM_LOCAL(recvbuf->buf, buf, ncopy);
    
    recvbuf->info[3] = false;
    
    nleft -= ncopy;
    buf   += ncopy;
    
    if (nleft) {
      buffer_number++;
      local_await(&recvbuf->info[3], buffer_number);
    }
  }
}



long MatchShmMessage(node, type)
     long node;
     long type;

{
   ShmemBuf *recvbuf;
   long  msg_type;

   recvbuf = TCGMSG_proc_info[node].recvbuf;

   if(recvbuf->info[3] == false) return (0); /* no message to receive */

   /* we have a message but let's see if want it */

   FLUSH_CACHE_LINE(recvbuf->info);
   COPY_FROM_LOCAL(recvbuf->info, &msg_type, sizeof(long));
   if(type == msg_type) return (1);
   return (0);
}
  
