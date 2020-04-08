

/*                                                    */
/*  $Author: mrdj $                                    */
/*  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $                        */
/*  $Locker:  $                                    */
/*  $Revision: 6317 $                                  */
/*  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/c.m,v $  */
/*  $State: Exp $                                     */
/*                                                    */

/*  ========== headers ============= */

/* compile SG with -ansi */
/* compile HP with -Ae */

#include <stdio.h>
#include <stdlib.h>      /* for exit, malloc, etc. */
#include <sys/stat.h>    /* for statcc */
#include <sys/types.h>   /* for lseek  */
#include <string.h>      /* for syserror() */
#include <errno.h>       /* use of perror */
#include <time.h>        /* for localtime */

/* this one generates warnings with pgcc */
#include <unistd.h>      /* various UNIX junk, POSIX, passwords, logins, etc */
#include <sys/times.h>  /* for etime */


#include <sys/time.h>    /* for gettimeofday */

#include <sys/fcntl.h>   /* for O_RDWR etc */


/*  ========== end of headers ============= */

/*  =========  local cpp constants ============= */

#include <errno.h>
#define	BLOCKSIZE	512*8	/* in bytes */
#define DEFMODE		0660	/* octal code for read/write group and owner */

#define OFF_T off_t
#define LSEEK lseek

/* extras */


/*  =========  Macros for fortran/C string passing ============= */

char *null_terminate(char *s, int len);
#define STRING_ARGS1 char *string1, INT *length1
#define STRING_ARGS2 char *string2, INT *length2
#define STRING_ARGS3 char *string3, INT *length3
#define GET_STRING1 (null_terminate(string1, *length1))
#define GET_STRING2 (null_terminate(string2, *length2))
#define GET_STRING3 (null_terminate(string3, *length3))


/*
		atmol io routines in c

	opencc ---- opens a file assigns a file descriptor
	putcc  ---- writes a 512 word block
	putccn  --- writes n 512 word blocks
	putccs  --- writes n 512 word blocks
	getcc  ---- reads a 512 word block
	getccn  --- reads n 512 word blocks
	getccs  ---- reads a 512 word block
	srchcc ---- positions the file at a given block
	srchccs ---- positions the file at a given block
	closecc --- closes files
        delcc   --- deletes (unlinks) the file
*/

/* ========== redefinitions for machines without  ========= */
/*              underscore append to fortran names          */


#ifdef CRAYXX
#define opencc_ OPENCC
#define srchcc_ SRCHCC
#define srchccs_ SRCHCCS
#define putcc_ PUTCC
#define putccn_ PUTCCN
#define putccs_ PUTCCS
#define getcc_ GETCC
#define getccn_ GETCCN
#define getccs_ GETCCS
#define closecc_ CLOSECC
#define delcc_ DELCC
#define iolog_ IOLOG
#endif

#ifdef NOC_
#define opencc_ opencc
#define srchcc_ srchcc
#define srchccs_ srchccs
#define putcc_ putcc
#define putccn_ putccn
#define putccs_ putccs
#define getcc_ getcc
#define getccn_ getccn
#define getccs_ getccs
#define closecc_ closecc
#define delcc_ delcc
#define iolog_ iolog
#endif


/* return type for c IO routines */
/* On MacOSX we have the following types:

   Type       #Bytes
   ----       ------
   int         4
   long        4
   long long   8

   So for an I8 build we need "long long".
   HvD, Aug 2009.

   */
#ifdef LONG_LONG_INTEGER
typedef  long long IORET;
typedef  long long INT;
typedef  long long UINT;
#elif LONG_INTEGER
typedef  long IORET;
typedef  long INT;
typedef  long UINT;
#else
typedef  int IORET;
typedef  int INT;
typedef  unsigned int UINT;
#endif

/*  =========  end local cpp constants ============= */

  
/* ===========  Function Prototypes/declarations ====================*/

void flushcc_(INT *unit,INT *ierror);
void opencc_(STRING_ARGS1, INT *unit,INT *ierror);
void srchcc_(INT *unit, INT *block, INT *ierror);
void srchccs_(INT *unit, INT *block, INT *ierror, INT *blsize);
void putcc_(INT *unit, char *buffer, INT *ierror);
void putccn_(INT *unit, char *buffer,INT *n,INT *ierror);
void putccs_(INT *unit, char *buffer,INT *ierror,INT *blsize);
void getcc_(INT *unit, char *buffer, INT *ierror);
void getccn_(INT *unit, char *buffer, INT *n, INT *ierror);
void getccs_(INT *unit, char *buffer, INT *ierror, INT *blsize);
void closecc_(INT *unit);
void delcc_(STRING_ARGS1, INT *ierror);
void syserr(char *msg);

/* ===========  end Function Prototypes ====================*/


/* file descriptors */

int fdss[110] = {0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0}; 

int mode = {DEFMODE};

/*    ======================  flushcc =====================  */
/*                  flush a file                             */

void flushcc_(INT *unit,        /* is pointer to the unit number */
	      INT *ierror)      /* error code */
{
/*      *ierror = fflush(fdss[*unit]); */
        *ierror = fflush(NULL);
	return;
}

/*    ======================  opencc =====================  */
/*                  opens a file for processing             */

void opencc_(STRING_ARGS1,       /* pointer to the file name  */
	     INT *unit,        /* is pointer to the unit number */
	     INT *ierror)      /* error code */
{
	static char path[200];
        *ierror = 0;
	strncpy(path,string1,*length1);
	path[*length1]='\0';
	if(( fdss[*unit] = open(path,(int)( O_RDWR|O_CREAT) , mode) ) == -1 ){
	  *ierror = 1;
          syserr("open"); 
	}
/*		fprintf(stderr,"%d %d\n",*unit,fdss[*unit]); */
	return;
}

void srchcc_(INT *unit,        /* is pointer to the unit number */
	     INT *block,       /* pointer to the block */
	     INT *ierror)      /* error code */
{
        *ierror = 0;

        if(LSEEK(fdss[*unit],(OFF_T)(*block-1)*BLOCKSIZE,0) == -1 ){
           *ierror = 1;
           syserr("search"); 
        }
	return;
}

void srchccs_(INT *unit,        /* is pointer to the unit number */
             INT *block,       /* pointer to the block */
             INT *ierror,      /* error code */
             INT *blsize)
{
        *ierror = 0;
/*      fprintf(stderr,"%d %d %d\n",*unit,fdss[*unit],*block); */
        if(lseek(fdss[*unit],(*block-1)*((*blsize)*8),0) == -1 ){
           *ierror = 1;
        }
        return;
}

/*  putcc: writes a record to a unit */
void putcc_(INT *unit,       /* pointer to the unit number */
	     char *buffer,    /* pointer to output from read */
	     INT *ierror)     /* error code */
{
	*ierror = 0;

    /* fprintf(stdout," about to output to unit %d fdss[unit]= %d BLOCKSIZE=%d \n", 
       *unit,fdss[*unit],BLOCKSIZE); */

        if(write(fdss[*unit],buffer,BLOCKSIZE) != BLOCKSIZE ){
	  *ierror = 1;
	  syserr("write");
	}
	return;
}

/*  putccn: writes a record to a unit */
void putccn_(INT *unit,      /* pointer to the unit number */
	    char *buffer,   /* pointer to output from read */
	    INT *n,         /* blocking factor */
	    INT *ierror)    /* error code */
{
  int nb;

  *ierror = 0;
  nb = (*n)*BLOCKSIZE;

  if(write(fdss[*unit],buffer,nb) != nb ){
    *ierror = 1;
    syserr("write");  
  }
  return;
}

void putccs_(INT *unit,       /* pointer to the unit number */
             char *buffer,    /* pointer to output from read */
             INT *ierror,     /* error code */
             INT *blsize)
{
        *ierror = 0;
        if(write(fdss[*unit],buffer,*blsize*8) != *blsize*8 ){
          *ierror = 1;
        }
        return;
}

void getcc_(INT *unit,     /* pointer to the unit number */
	    char *buffer,  /* is pointer to output from read */
	    INT *ierror)   /* error code */
{
  *ierror = 0;
  if(read(fdss[*unit],buffer,BLOCKSIZE) != BLOCKSIZE ){
    *ierror = 1;  
    syserr("read"); 
   }
  return;
}
/* getccn: reads a record form a unit */
void getccn_(INT *unit,    /* pointer to the unit number */
	    char *buffer,   /* is pointer to output from read */
	    INT *n,         /* is the number of blocks to  be read */
	    INT *ierror)    /* error code */
{
  int nb;

  *ierror = 0;
  nb = (*n)*BLOCKSIZE;   

  if (read(fdss[*unit],buffer,nb) != nb ){
    *ierror = 1;
    syserr("read");  
  }
  return;
}
/* getccs: reads a record form a unit */
void getccs_(INT *unit,     /* pointer to the unit number */
            char *buffer,  /* is pointer to output from read */
            INT *ierror,   /* error code */
            INT *blsize)
{
  *ierror = 0;
  if(read(fdss[*unit],buffer,*blsize*8) != *blsize*8 ){
    *ierror = 1;
   }
  return;
}

/* ================ closcc: closes a unit =============  */

void closecc_(INT *unit)
{
  if(fdss[*unit] != 0) close(fdss[*unit]);
  return;
}

/*  ===============  delcc: deletes a file =============  */

void delcc_(STRING_ARGS1,    /*  pointer to the file name */
	     INT *ierror)      /* error code */
{
        static char path[200];
	*ierror = 0;
	strncpy(path,string1,*length1);
	path[*length1]='\0';
	if(unlink(path) < 0 ){
	  *ierror = 1;
          printf("problem with deletion -%s- \n",path);
	  syserr("delete"); 
	}
	return;
}

/*  ============  syserr: print the error message ====== */
/* (internal to this file and not called from FORTRAN */

void syserr(char *msg)
{
/*  extern int errno, sysnerr; */

/* 
   code for UNIX machines where we can easily put the
   character string back into fortran 
*/
  char msgbuf[101];

  if(errno > 0 )sprintf(msgbuf,"%-29s", strerror(errno));
  else sprintf(msgbuf,"%-49s","(no system code)");
  (void) iolog_(&errno, msgbuf);
/*	exit(1);  */
}

/* C packing */

void pck8c_(long *ints, unsigned char *chars, long *pn);
void upck8c_(long *ints, unsigned char *chars, long *pn);

#ifdef NOC_
#define  pck8c_ pck8c
#define  upck8c_ upck8c
#endif

#ifdef CRAYXX
#define  pck8c_ PCK8C
#define  upck8c_ UPCK8C
#endif

/* ========== pck8c: copies n long integer into character ======== */
/* modification for convex fortran interface  _copy */

void pck8c_(long *ints, unsigned char *chars, long *pn)
{
  long n,*ints_copy;
  int b;
  unsigned char *chars_copy;
  n= *pn;
  ints_copy = ints;
  chars_copy = chars;
  if (n==0) return;
  for (b=1;b<=n;b++) *(chars_copy++) = *(ints_copy++);
  return;
}

/* ========= upck8c: copies n characters into integers ========== */
/* modification for convex fortran interface  _copy */

void upck8c_(long *ints, unsigned char *chars, long *pn)
{
  long n,*ints_copy;
  int b;
  unsigned char *chars_copy;
  n= *pn;
  ints_copy = ints;
  chars_copy = chars;
  if (n==0) return;
  for (b=1;b<=n;b++) {
    *(ints_copy++) = *(chars_copy++);
  }
  return;
}

/* ================ Utility routines ===================  */

#ifdef NOC_
#define cordmp_ cordmp
#define statcc_ statcc
#define daycc_ daycc
#define namcc_ namcc
#define hostcc_ hostcc
#define pidcc_ pidcc
#define uidcc_ uidcc
#define gtnvcc_ gtnvcc
#define walltimecc_ walltimecc
#define cputimecc_ cputimecc
#define ver_c_ ver_c
#define mallocc_ mallocc
#define gmem_c_pointer_diff_ gmem_c_pointer_diff
#define mallocc2_ mallocc2
#define freec_  freec
#define chkadr_  chkadr
#endif

#ifdef LINUXF2C
#define ver_c_ ver_c__ 
#define gmem_c_pointer_diff_ gmem_c_pointer_diff__
#endif

#ifdef GFS
#define ver_c_ ver_c_
#define gmem_c_pointer_diff_ gmem_c_pointer_diff_
#endif

#ifdef CRAYXX
#define cordmp_ CORDMP
#define statcc_ STATCC
#define daycc_ DAYCC
#define namcc_ NAMCC
#define walltimecc_ WALLTIMECC
#define cputimecc_ CPUTIMECC
#define hostcc_ HOSTCC
#define pidcc_ PIDCC
#define uidcc_ UIDCC
#define gtnvcc_ GTNVCC
#define ver_c_ VER_C
#define mallocc_ MALLOCC
#define gmem_c_pointer_diff_ GMEM_C_POINTER_DIFF
#define mallocc2_ MALLOCC2
#define freec_  FREEC
#define chkadr_  CHKADR
#endif

/* ===== cordmp: core dump ========== */
void cordmp_()
{
    abort();
}

/* ======= statcc : File status/size =========== */
static struct stat buff; 
void statcc_(STRING_ARGS1,
	   INT *isize)
{
   int iret;
   char *file;
   file = GET_STRING1;
   iret = stat(file,&buff);
   free(file);
   if(iret == -1){
     *isize=-1;
   }else{
     *isize = (int) buff.st_size;
   }
}

/* ======= daycc : Date string =========== */
void daycc_(STRING_ARGS1)
{

#if defined _SYSTYPE_SVR4
   int tzp;
#else
#if ! defined NEC
  struct timezone tzp;
#endif
#endif
  struct timeval tp;
  char *ctime();
  time_t i;
#if defined NEC
  gettimeofday(&tp);
#else
  gettimeofday(&tp,&tzp);
#endif
  i = tp.tv_sec;
  strncpy(string1,ctime(&i),24);
  return;
}

/* ========  namcc : UNIX Username  ========== */
void namcc_(STRING_ARGS1)
{
  char buffer[32];
  char *t, *cuserid();
  int l, i;

  t = cuserid(NULL);
  if(t){
    strncpy(buffer,t,31);
  }else{
    strcpy(buffer,"n/a");
  }
  l=strlen(buffer);
  if(l > *length1){
    strncpy(string1,buffer,*length1);
  }else{
    strcpy(string1,buffer);
    for(i=l;i<*length1;i++)*(string1+i)=' ';
  }
  return;
}

/**************************************************************************

                                 t i m e r s

 timers for wall-clock  (walltimecc)  and cpu time (cputimecc) 

 return real*8 values in seconds 

 real*8 tbuf(3)
 call cputimecc(tbuf)
 c  tbuf(1) - user cpu
 c  tbuf(2) - system cpu
 c  tbuf(3) - total cpu

 In some cases a fortran-called system routine is used instead (see machscf.m)

 **************************************************************************/

static int iflag=0;
void walltimecc_(double *elp)
{
  static struct timeval tp, tp0;
#if defined _SYSTYPE_SVR4
  int tzp;
#else
#if ! defined NEC
  struct timezone tzp;
#endif
#endif
  double t;
  if(iflag){
#if defined NEC
    gettimeofday(&tp);
#else
    gettimeofday(&tp,&tzp);
#endif
    t =  (double) (tp.tv_sec - tp0.tv_sec);
    *elp = t;
    t = (double) (tp.tv_usec - tp0.tv_usec);
    *elp += (double) t / 1000000.0 ;
  }else{
#if defined NEC
    gettimeofday(&tp0);
#else
    gettimeofday(&tp0,&tzp);
#endif
    iflag=1;
    *elp=0.0;
  }
}

void cputimecc_(double *cpu)
{
   struct tms timestruct;
   long iret, sec, sec1;
   static long tick = 0;
   if(!tick)tick = sysconf(_SC_CLK_TCK);
   iret = times(&timestruct);
   sec = timestruct.tms_utime/tick;
   sec1 =  timestruct.tms_utime - sec*tick;
   cpu[0] = (double) sec + (double) sec1 / (double) tick;
   sec = timestruct.tms_stime/tick;
   sec1 =  timestruct.tms_stime - sec*tick;
   cpu[1] = (double) sec + (double) sec1 / (double) tick;
   cpu[2] = cpu[0] + cpu[1];
}
/**************************************************************************

    Calls for system information

 **************************************************************************/

/* ========  hostcc : UNIX Hostname  ========== */
void hostcc_(STRING_ARGS1)
{
  char buffer[32];
  int l, i;
  gethostname(buffer,32);
  l=strlen(buffer);
  if(l > *length1){
    strncpy(string1,buffer,*length1);
  }else{
    strcpy(string1,buffer);
    for(i=l;i<*length1;i++){
      *(string1+i)=' ';
    }
  }
  return;
}
/* ========  pidcc : UNIX Process id  ========== */
void pidcc_(INT *i)
{
  *i = (INT) getpid();
  return;
}
/* ========  uidcc : UNIX User id  ========== */
void uidcc_(INT *i)
{
  *i = (INT) getuid();
  return;
}



/********************************************************************
 * null_terminate(s,n) null terminate a string from FORTRAN, and copy it
 ********************************************************************/
char *null_terminate(char *s, int n)
{
  int i;
  char *t;
  for(i = n - 1; *(s+i) == ' ' && i; i--)
    ;
  t = (char *) malloc(sizeof(char) * (i+2));
  if(!t){
     fprintf(stderr,"Memory allocation error\n");
     abort();
     exit(1);
  }
  strncpy(t,s,i+1);
  *(t+i+1) = '\0';
  return t;
}

/* ================== getmem :  memory allocation ==============

   getmem gets n real*8 storage locations and returns its
   address (iaddr) and offset (ioff) within the real*8 array work
   so that the usable memory is (work(i+ioff),i=1,n).
   e.g. 
        call getmem(n,work,iaddr,ioff)
        if (iaddr.le.0 .or. ioff.le.0) call error

   Mods are needed to release this later. */

#ifdef NOC_
#define getmem_ getmem
#endif
#ifdef CRAYXX
#define getmem_ GETMEM
#endif
void getmem_(INT *pn, double *pwork, UINT *paddr, UINT *pioff);

void getmem_(INT *pn, 
	    double *pwork, 
	    UINT *paddr, 
	    UINT *pioff)
{
 unsigned int size = 8;
 int iret;
 
 union screwup {
   unsigned long integer;
   double *address;
 } fiddle;
 unsigned long offset;

 /*iret = posix_memalign(&ptemp, size, (unsigned) size* *pn); */
 /* printf("checking sizes %d %d\n", sizeof(unsigned long), sizeof (double *)); */

 if (*pn < 0) {
    /* Oh shit, something has suffered an overflow in the Fortran. 
       Bail out now!!! */
    *paddr = 0;
    *pioff = 0;
    return;
 }
 fiddle.address = (double *) malloc( size * *pn );
/* printf("allocating %ld\n",fiddle.address); */
 offset = fiddle.integer & 7;

 if (offset != 0) {
   printf("Aligment error in getmem\n");
   fprintf(stderr,"Aligment error in getmem\n");
   exit(-1);
 }
 *paddr = (UINT) fiddle.address;
 *pioff = (UINT) (fiddle.address - pwork);
}


/* prpnam: converts string like cnv and number 5 to string c5v */
/* strings must not end with \0, because fortran .eq. would not work */
/* simplified by PS for ipsc where isupper/tolower failed */

#ifdef NOC_
#define prpnam_ prpnam
#endif

#ifdef CRAYXX
#define prpnam_ PRPNAM 
#endif

void prpnam_(STRING_ARGS1, STRING_ARGS2, INT *order)
{
  register int i;
  strncpy(string1,string2,*length1);
  if(!strncmp(string1,"s2n",3) || !strncmp(string1,"S2N",3) ){
    string1[1] = '0'+2*(*order);
    string1[2] = ' ';
  }else{
    for(i=0;i<*length1;i++){
      if(string1[i]=='n' || (string1[i]=='N'))string1[i]= '0'+(*order); 
    }
  }
}
void gtnvcc_(STRING_ARGS1,STRING_ARGS2,INT *ier)
{
  char *getenv(), *null_terminate(), *p, *p1;
  int i,l;
  p1 = GET_STRING1;
  p = getenv(p1);
  free(p1);
  if(!p){
     *ier = -1;
     for(i=0;i<*length2;i++)*(string2+i)=' ';
     return;
  }
  l = strlen(p);
  if(l > *length2 ){
     strncpy(string2,p,*length2);
     *ier = 1;
  }else{
     strcpy(string2,p);    
     for(i=l;i<*length2;i++)*(string2+i)=' ';
     *ier = 0;
  }
}

/*  ==== end of ipsc special routines ===== */



/* used in gmem.M */

void gmem_c_pointer_diff_(double *q1, double *q2, INT *offset)
{
  *offset = (INT)  (q1 - q2);
}

static char source[]="$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/c.m,v $";
static char revision[]="$Revision: 6317 $";
static char date[]="$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $";

void ver_c_(STRING_ARGS1,STRING_ARGS2,STRING_ARGS3)
{
  int l, start;
  l = strlen(source);
  start = 8;
  l-=start;
  if(l > *length1)l = *length1;
  strncpy(string1,source+start,l);

  l = strlen(revision);
  start = 10;
  l-=start;
  if(l > *length2)l = *length2;
  strncpy(string2,revision+start,l);

  l = strlen(date);
  start = 6;
  l-=start;
  if(l > *length3)l = *length3;
  strncpy(string3,date+start,l);
}  


/* mallocc */

void mallocc_(void *q1, void *q2, INT *size, INT *off, INT *handle)
{
  void *p;
  long offset, realsize;

  realsize =  (long) q2 - (long) q1;

/*  printf("request %ld objects of size %ld bytes\n",*size, realsize);*/

  p = (void *) malloc( *size * realsize );

/*  printf("allocated -address %ld \n", p );*/
   
  if(!p){
    *off = 0;
    return;
  }
  offset = (long ) p - (long) q1;

/*  printf("allocated -address %ld  offset %ld bytes\n", p, offset);*/
  fflush(stdout);


  *off = offset / realsize;
  *handle = (INT) ( (long) p - (long) q1 );

}

void mallocc2_(void *q1, void *q2, INT *size, INT *off, INT *handle)
{
  void *p;
  long offset, realsize;

  realsize =  (long) q2 - (long) q1;

  p = (void *) malloc( *size * realsize );

  if(!p){
    *off = 0;
    return;
  }

  offset = (long ) p - (long) q1;

  /*printf("malloc: allocate %ld \n", p); */

  *handle = (INT) ( (long) p - (long) q1 );


  *off = 1 + offset / realsize;  /* add 1 so that qq(ivoff+*ioff) means what we want */

  /*  printf("check allignment.. \n p=%ld \n q1=%ld \n offset=%ld \n *off = %ld \n q1+*off*8= %ld\n",p,q1,offset,*off, (long) q1 + 8 * *off);*/

}
void freec_(void *q1, INT *handle, INT *stat)
{

  /*  printf("malloc: free %ld \n", ((long) q1 + (long) *handle));*/

  free((void *) ((long) q1 + (long) *handle));

  *stat = 0;

}
void chkadr_(void *q1)
{
  /* printf("malloc: chkadr %ld \n", q1);*/
}

void chkadr2_(void *q1, INT *addr)
{
  *addr = (INT) q1;
}



/*

   C-implementations of these routines use unsigned integers and 
   so allow larger numbers to be packed.

*/
#ifdef CRAYXX
#define pack3cc_ PACK3CC
#define upack3cc_ UPACK3CC
#endif

#ifdef NOC_
#define pack3cc_ pack3cc
#define upack3cc_ upack3cc
#endif

void pack3cc_(unsigned short *stuff,
	INT *i,
	INT *j,
	INT *k){

  stuff[0] = *i/65536;
  stuff[1] = *i-stuff[0]*65536;
  stuff[2] = *j;
  stuff[3] = *k;
}

void upack3cc_(unsigned short *stuff,
	  INT *i,
	  INT *j,
	  INT *k){

  *i = stuff[0]*65536+stuff[1];
  *j = stuff[2];
  *k = stuff[3];
}



#ifdef CRAYXX
#define exitc_ EXITC
#endif

#ifdef NOC_
#define exitc_ exitc
#endif

void exitc_(code)
     INT *code;
{
  (void) fflush(stdout);
  (void) fflush(stderr);
  exit(*code);
}


#ifdef CRAYXX
#define stderrc_ STDERRC
#endif

#ifdef NOC_
#define stderrc_ stderrc
#endif

void stderrc_(STRING_ARGS1,
	     INT *code,
	     INT *sys,
	     INT *node)
{
  char *t;
  t = GET_STRING1;

  (void) fflush(stdout);
  (void) fflush(stderr);

  if(*code)(void) fprintf(stderr, "%3d: GAMESS-UK Error %ld : %s \n", *node, *code, t);
  else(void) fprintf(stderr, "%3d: GAMESS-UK Error: %s \n", *node, t);

  if(*sys) (void) perror("system message");

  free(t);

  (void) fflush(stdout);
  (void) fflush(stderr);
}

#ifdef CRAYXX
#define ishftcc_ ISHFTCC
#endif

void ishftcc_(UINT *in, 
	      UINT *out,
    	      INT *bits)
{
  if(*bits  < 0){
    *out = *in>>-(*bits);
  }else{
    *out = *in<<*bits;
  }
}
void btestcc_(UINT *in, 
	      UINT *out,
	      INT *bit)
{
  UINT mask;
  mask=1<<(*bit);
  *out = (INT) (*in &mask);
}
void ibsetcc_(UINT *in, 
	 UINT *out,
	 INT *bit)
{
  UINT mask;
  mask=1<<(*bit);
  *out = *in | mask;
}
void ibclrcc_(UINT *in, 
	 UINT *out,
	 INT *bit)
{
  UINT mask;
  mask = 0xfffffffe<<(*bit);
  *out = *in & mask;
}


/* Added for alex turners matrix code */
void flushout_() {
  fflush(stdout);
}



void segfaulthandler(int sig) {
void *trace[100];
exit(1);
}

void siginthandler(int sig) {
void *trace[100];
exit(1);
}

/* After this segfaulthandler will be called */
void installsegfaulthandler_() {
}

/* After this the default segv handler will be called */
void uninstallsegfaulthandler_() {
}

/* Create a stack trace by Interrupting (ctrl-c) a hanging job */
void installctrlchandler_() {
}

/* reinstate default ctrl-c */
void uninstallctrlchandler_() {
}

void fgetdate_(char zdate[])
{
    time_t timet;

    timet=time(NULL);
    strncpy(zdate,ctime(&timet),24);

}
