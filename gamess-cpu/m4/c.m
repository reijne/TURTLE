_IF(win95)
_IFN(ifc)
#pragma aux srchcc "^"
#pragma aux srchccs "^"
#pragma aux putcc "^"
#pragma aux putccn "^"
#pragma aux putccs "^"
#pragma aux getcc "^"
#pragma aux getccn "^"
#pragma aux getccs "^"
#pragma aux closecc "^"
#pragma aux timef "^"
#pragma aux pck8c "^"
#pragma aux upck8c "^"
#pragma aux cordmp "^"
#pragma aux getmem "^"
#pragma aux memalign "^"
#pragma aux cputm "^"
#pragma aux setcw "^"
#pragma aux sett0tijd "^"
#pragma aux upack3cc "^"
#pragma aux pack3cc "^"
#pragma aux ver_c "^"
#pragma aux ver_tsortc "^"
#pragma aux walltimecc "^"
#pragma aux cputimecc "^"
#pragma aux mallocc "^"
#pragma aux gmem_c_pointer_diff "^"
#pragma aux daycc "^"
#pragma aux namcc "^"
#pragma aux uidcc "^"
#pragma aux pidcc "^"
#pragma aux ishftcc "^"
#pragma aux ibsetcc "^"
#pragma aux ibclrcc "^"
#pragma aux btestcc "^"
#pragma aux freec "^"
#pragma aux setcw "^"
#define NOC_
#define ishftcc_ ishftcc
#define ibsetcc_ ibsetcc
#define ibclrcc_ ibclrcc
#define btestcc_ btestcc
#define gtnvcc_ gtnvcc
_ENDIF
_ENDIF

_IF(taskfarm)
#include <setjmp.h>
_ENDIF

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

_IFN(win95)
/* this one generates warnings with pgcc */
#include <unistd.h>      /* various UNIX junk, POSIX, passwords, logins, etc */
#include <sys/times.h>  /* for etime */
_ENDIF

_IF(win95)
#include <io.h>          /* for lseek, open, close, etc. */
#include <process.h>     /* for getpid */
_ENDIF

_IF(win95,hp700)
#include <time.h>        /* for gettimeofday */
#include <sys/timeb.h>
_ELSE
#include <sys/time.h>    /* for gettimeofday */
_ENDIF

_IF(rs6000,t3d,unicos,win95,hp700)
#include <fcntl.h>       /* for O_RDWR etc */
_ELSE
#include <sys/fcntl.h>   /* for O_RDWR etc */
_ENDIF

_IF(glibc)
#include <execinfo.h>    /* for backtrace, backtrace_symbols_fd */
_ENDIF

/*  ========== end of headers ============= */

/*  =========  local cpp constants ============= */

#include <errno.h>
#define	BLOCKSIZE	512*8	/* in bytes */
_IFN(vms)
#define DEFMODE		0660	/* octal code for read/write group and owner */
_ELSE
#define DEFMODE		0660	/* octal code for read/write group and owner */
_ENDIF

_IF(mips4)
#define OFF_T off64_t
#define LSEEK lseek64
_ELSE
#define OFF_T off_t
#define LSEEK lseek
_ENDIF

/* extras */
_IF(mips4)
off64_t lseek64();
_ENDIF

_IF(ksr)
#define int long 
_ENDIF

/*  =========  Macros for fortran/C string passing ============= */

_IF(t3d)
char *null_terminate(char *s, int len);
#define STRING_ARGS1 char *string1, int string1_length, int *length1
#define STRING_ARGS2 char *string2, int string2_length, int *length2
#define STRING_ARGS3 char *string3, int string3_length, int *length3
_ELSE
char *null_terminate(char *s, int len);
#define STRING_ARGS1 char *string1, INT *length1
#define STRING_ARGS2 char *string2, INT *length2
#define STRING_ARGS3 char *string3, INT *length3
_ENDIF
#define GET_STRING1 (null_terminate(string1, *length1))
#define GET_STRING2 (null_terminate(string2, *length2))
#define GET_STRING3 (null_terminate(string3, *length3))


_IF(cio)
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
_IF(win95)
_IF(ifc)
#define gmem_c_pointer_diff_ GMEM_C_POINTER_DIFF
#define prpnam_ PRPNAM
#define ver_c_ VER_C
#define pack3cc_ PACK3CC
#define upack3cc_ UPACK3CC
#define flushcc_ FLUSHCC
#define flushout_ FLUSHOUT
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
#define ishftcc_ ISHFTCC
#define ibsetcc_ IBSETCC
#define ibclrcc_ IBCLRCC
#define btestcc_ BTESTCC
#define gtnvcc_ GTNVCC
#define exitc_ EXITC
#define walltimecc_ WALLTIMECC
#define cputimecc_ CPUTIMECC
#define pidcc_ PIDCC
#define statcc_ STATCC
#define stderrc_ STDERRC
#define cordmp_ CORDMP
#define getmem_ GETMEM
#define mallocc_ MALLOCC
#define freec_ FREEC
#define setcw SETCW
#define getdate GETDATE
_ENDIF
_ENDIF

_IF(rs6000)
_IF(rs6000_noextname)
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
#define hostcc_ hostcc
#define namcc_ namcc
#define exitc_ exitc
#define pidcc_ pidcc
#define gmem_c_pointer_diff_ gmem_c_pointer_diff
#define statcc_ statcc
#define getmem_ getmem
#define upack3cc_ upack3cc
#define pack3cc_ pack3cc
#define ver_c_ ver_c
#define mallocc_ mallocc
#define freec_ freec
#define daycc_ daycc
#define uidcc_ uidcc
#define gtnvcc_ gtnvcc
#define stderrc_ stderrc
#define cordmp_ cordmp
#define etime_ etime
#define prpnam_ prpnam
_ENDIF(rs6000_noextname)
_ENDIF(rs6000)

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
_IF(alliant)
typedef  void  IORET;
typedef  int INT;
_ELSE
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
_ENDIF

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
_IFN(vms)
_IF(ksr)
	if(( fdss[*unit] = open(path,(long)( O_RDWR|O_CREAT) , mode) ) == -1 ){
_ELSE
_IF(win95)
	if(( fdss[*unit] = open(path,(int)( O_BINARY|O_RDWR|O_CREAT) , mode) ) == -1 ){
_ELSE
	if(( fdss[*unit] = open(path,(int)( O_RDWR|O_CREAT) , mode) ) == -1 ){
_ENDIF
_ENDIF
_ELSE
        /* fprintf(stderr,"filename length is %d \n",*length1); */
        /* fprintf(stderr,"opencc for %d %d %s\n",*unit,fdss[*unit],path); */
	if(( fdss[*unit] = open(path,(int)( O_RDWR|O_CREAT) , mode,"shr=upi",
            "fop=cif") ) == -1 ){
_ENDIF
	  *ierror = 1;
_IF(vms)
/*        if(( fdss[*unit] = open(path,(int)( O_RDWR|O_CREAT),mode) ) == -1){
          ierror = 1;*/		  
	  perror("open");
_ENDIF
          syserr("open"); 
	}
_IFN(vms)
/*		fprintf(stderr,"%d %d\n",*unit,fdss[*unit]); */
_ELSE
/*        fprintf(stderr,"%d %d\n",*unit,fdss[*unit]);*/ 
/*          printf("opencc for %d %d %s\n",*unit,fdss[*unit],path);*/
_ENDIF
	return;
}

void srchcc_(INT *unit,        /* is pointer to the unit number */
	     INT *block,       /* pointer to the block */
	     INT *ierror)      /* error code */
{
        *ierror = 0;

_IF(apollo,convex,alliant)
_IF(vms)
	if(lseek(fdss[*unit],(*block-1)*BLOCKSIZE,SEEK_SET) == -1 ){
_ELSE
	if(lseek(fdss[*unit],(*block-1)*BLOCKSIZE,L_SET) == -1 ){
_ENDIF
_ELSE
        if(LSEEK(fdss[*unit],(OFF_T)(*block-1)*BLOCKSIZE,0) == -1 ){
_ENDIF
           *ierror = 1;
_IF(vms)
	   perror("search");
_ENDIF
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
_IF(vms)
         int iwritten;
_ENDIF
	*ierror = 0;
_IFN(vms)

    /* fprintf(stdout," about to output to unit %d fdss[unit]= %d BLOCKSIZE=%d \n", 
       *unit,fdss[*unit],BLOCKSIZE); */

        if(write(fdss[*unit],buffer,BLOCKSIZE) != BLOCKSIZE ){
	  *ierror = 1;
	  syserr("write");
_ELSE

/*  fprintf(stdout," about to output to unit %d fdss[unit]= %d BLOCKSIZE=%d \n",
	  *unit,fdss[*unit],BLOCKSIZE); */

/*	if(write(fdss[*unit],buffer,BLOCKSIZE) != BLOCKSIZE ){*/
	iwritten = write(fdss[*unit],buffer,BLOCKSIZE);
        if (iwritten != BLOCKSIZE){
        fprintf(stderr," %d bytes written from write \n",iwritten);
        *ierror = 1;
	  perror("write");
          syserr("write");
_ENDIF
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
_IFN(vms)
    syserr("write");  
_ELSE
    perror("write");
    syserr("write");  
_ENDIF
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
_IFN(vms)
    syserr("read"); 
_ELSE
    perror("read");
    syserr("read"); 
_ENDIF
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
_IFN(vms)
    syserr("read");  
_ELSE
    perror("read");
    syserr("read");  
_ENDIF
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
_IFN(vms)
	if(unlink(path) < 0 ){
_ELSE
	if(delete(path) < 0 ){
_ENDIF
	  *ierror = 1;
_IF(vms)
	  perror("delete");
_ENDIF
          printf("problem with deletion -%s- \n",path);
	  syserr("delete"); 
	}
	return;
}

/*  ============  syserr: print the error message ====== */
/* (internal to this file and not called from FORTRAN */

void syserr(char *msg)
{
_IF(ksr)
  extern long sysnerr;
_ELSEIFN(vms)
/*  extern int errno, sysnerr; */
_ELSE
  extern int errno;
/*, sysnerr;*/
_ENDIF

_IF(unix)
/* 
   code for UNIX machines where we can easily put the
   character string back into fortran 
*/
  char msgbuf[101];

  if(errno > 0 )sprintf(msgbuf,"%-29s", strerror(errno));
  else sprintf(msgbuf,"%-49s","(no system code)");
  (void) iolog_(&errno, msgbuf);
_ELSE
  /* write message to stderr */
	fprintf(stderr,"Error in atmol i/o %s (%d ",msg,errno);
_IFN(vms)
	if(errno > 0 )
		fprintf(stderr,"; %s )\n", strerror(errno));
_ELSE
	if(errno > 0 ){
        /* fprintf(stderr,"; %s )\n", strerror(errno));}*/
	  fprintf(stderr," error in i/o ");}
_ENDIF
	else
		fprintf(stderr,")\n");
_ENDIF
/*	exit(1);  */
_IF(vms)
  return (0);
_ENDIF
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
_IF(ksr)
  long b;
_ELSE
  int b;
_ENDIF
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
_IF(ksr)
  long b;
_ELSE
  int b;
_ENDIF
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
_ELSE
/* repeat typedef for non cio */
_IF(alliant)
typedef  void  IORET;
typedef  int INT;
_ELSE
#ifdef LONG_INTEGER
typedef  long IORET;
typedef  long INT;
typedef  long UINT;
#else
typedef  int IORET;
typedef  int INT;
typedef  unsigned int UINT;
#endif
_ENDIF

_ENDIF(cio)

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

_IF(rs6000,hp700,hpux11,t3d,unix)
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
_ENDIF

_IFN(win95)
_IF(rs6000,hp700,hpux11,t3d,nec,unix)
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
_ENDIF

_IF(rs6000,hp700,hpux11,parallel,unix)
_IFN(macosx)
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
_ENDIF
_ENDIF

_IFN(ipsc,t3d)
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
_ENDIF

_IF(t3d,ipsc)
/* cputimecc not suitable for t3d due to 4-byte float - use tsecnd */
/*           on ipsc, use dclock */
_ELSE
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
_ENDIF
/**************************************************************************

    Calls for system information

 **************************************************************************/

_IF(rs6000,hp700,hpux11,parallel,unix,itanium,opteron,xeon,em64t)
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
_ENDIF

_IF(cray,unicos)
/* ======== linkcc:  create a symbolic link (UNICOS) ====== */ 
/* #include <stdio.h> */
LINKCC(iunit,link,len)
     int *iunit;
     char *link;
     int *len;
{
  int ii;
  char *l1, *null_terminate();
  char fname[30];
  sprintf(fname,"fort.%d",*iunit);
  l1 = null_terminate(link,*len);
  if(ii = symlink(l1,fname)){
      fprintf(stderr,"symlink failed %d\n",ii);
  }   
  free(l1);
  return 0;
}
_ENDIF

_ELSE
/* Windows 95 timers ( before this till ifn(win95) skipped */

void pidcc_(INT *i)
{
  *i = (INT) getpid();
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
  struct timeb t;
  static struct timeb t0;

  if(iflag){
    ftime(&t);
    *elp=(double)(t.time-t0.time)+(double)(t.millitm-t0.millitm)/1000;
  }else{
    ftime(&t0);
    iflag=1;
    *elp=0.0;
  }
}

void cputimecc_(double *cpu)
{
  clock_t t;
  static clock_t t0;
  static int iflagj=0;

  if(iflagj){
    t=clock();
    *cpu=(double)(t-t0)/CLOCKS_PER_SEC;
  }else{
    t0=clock();
    iflagj=1;
    *cpu=0.0;
  }
}

void getdate(char zdate[])
{
    struct tm tp;
    time_t timet;

    timet=time(NULL);
_IF(ifc)
    tp = *localtime(&timet);
_ELSE
    _localtime(&timet,&tp);
_ENDIF
    strncpy(zdate,asctime(&tp),24);

}

#include <float.h>

#define IEEE  (_IC_AFFINE | _RC_NEAR | _PC_64 \
               |_EM_INVALID |_EM_DENORMAL|_EM_OVERFLOW\
               |_EM_ZERODIVIDE|_EM_UNDERFLOW)
                      
  void setcw()
  { unsigned int fp_cw = 0 ;
  
    fp_cw =_control87(IEEE,IEEE);
  }

_ENDIF

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

_IF(convex)
/* #include <sys/time.h> */
/* #include <sys/resource.h> */
#define	RUSAGE_SELF	0
#define RUSAGE_CHILDREN	(-1)

void cgetru_(who,rusage)
INT  who, *rusage ;
{
	getrusage(who,rusage);
	return;
}
_ENDIF
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
_IF(64bitpointers)
void getmem_(INT *pn, double *pwork, UINT *paddr, long *pioff);

void getmem_(INT *pn,
            double *pwork,
            UINT *paddr,
            long *pioff)
_ELSE
void getmem_(INT *pn, double *pwork, UINT *paddr, UINT *pioff);

void getmem_(INT *pn, 
	    double *pwork, 
	    UINT *paddr, 
	    UINT *pioff)
_ENDIF
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
_IF(64bitpointers)
 *pioff = (long) (fiddle.address - pwork);
_ELSE
 *pioff = (UINT) (fiddle.address - pwork);
_ENDIF
}
_IF(cray,unicos)
/* #include <sys/times.h> */
#include <sys/machd.h> 
void USTIME(float *out)
/* Return User+System CPU time (in secs) of process and its children 
   and use instead of second to be consistent with nqs
   see man times for INFO
   Andries de Man, Utrecht 1992 */
{ struct tms buffer;
  times(&buffer);
  *out = ((float) buffer.tms_utime+buffer.tms_stime+
                        buffer.tms_cutime+buffer.tms_cstime)/HZ;
  return;
}
_ENDIF 

_IF(vms)
#include <types.h> 
#include <socket.h> 
_ENDIF

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
_IF(ipsc,t3d,unix)
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
_ENDIF

/*  ==== end of ipsc special routines ===== */

_IF(sun)
/*  =========== catchit:  signal handler for sun */
#include <floatingpoint.h> 
#include <signal.h>
static void catchit()
{
  printf("!!  Floating point interrupt caught  !!\n");
  fflush(stdout);
  (void) signal(SIGIOT, SIG_DFL);
  abort();
}
void ieeetrap_()
{
/* (void) ieee_handler("set","inexact",catchit);*/
  (void) ieee_handler("set","underflow",SIGFPE_ABORT);
}
_ENDIF

_IF(dec)
_IF(vms)

int fgethostname(char *host, INT len);
INT fgetenv(char *ename, char *evalue);

INT fgethostname(char *host, INT len)
{
   return(gethostname(host, len));
}
INT fgetenv(char *ename, char *evalue)
{
   evalue=(getenv(ename));
   return(1);
}

_ENDIF
_ENDIF

/* used in gmem.M */

_IF(64bitpointers)
void gmem_c_pointer_diff_(double *q1, double *q2, long *offset)
{
  *offset = (long)  (q1 - q2);
}
_ELSE
void gmem_c_pointer_diff_(double *q1, double *q2, INT *offset)
{
  *offset = (INT)  (q1 - q2);
}
_ENDIF

_IF(t3d)
FILE *fplog;
INITLOG(int *proc)
{
   char name[20];
   sprintf(name,"trace.%d",*proc);
   fplog =fopen(name,"w");
}
LOGSEND(int *dest, int *src, int *len, int *tag,STRING_ARGS1)
{
   char *key;
   key=GET_STRING1;
   fprintf(fplog,"%4d %4d %4d %4d send %s\n",*dest,*src, *len, *tag, key);
   free(key);
   fflush(fplog);
}
LOGRCV(int *dest, int *src, int *len, int *tag,STRING_ARGS1)
{
   char *key;
   key=GET_STRING1;
   fprintf(fplog,"%4d %4d %4d %4d rcv %s\n",*dest,*src, *len, *tag, key);
   free(key);
   fflush(fplog);
}
LOGSYNC(int *tag,STRING_ARGS1)
{
   char *key;
   key=GET_STRING1;
   fprintf(fplog,"            %4d synch %s\n",*tag, key);
   free(key);
   fflush(fplog);
}
_ENDIF
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

_IF(unix,t3d,mips4,cray,dec)

/* mallocc */

_IF(64bitpointers)
void mallocc_(void *q1, void *q2, INT *size, long *off, long *handle)
_ELSE
void mallocc_(void *q1, void *q2, INT *size, INT *off, INT *handle)
_ENDIF
{
  void *p;
  long offset, realsize;

_IF(cray)
  /* for true cray systems word length computes to 1, hardwire to 8 */
  realsize = 8;
_ELSE
  realsize =  (long) q2 - (long) q1;
_ENDIF

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

_IF(cray)
 /* cray offset computes directly in words */
  realsize=1;
_ENDIF

_IF(64bitpointers)
  *off = offset / realsize;
  *handle = (long) p - (long) q1;
_ELSE
  *off = offset / realsize;
  *handle = (INT) ( (long) p - (long) q1 );

_ENDIF
}

void mallocc2_(void *q1, void *q2, INT *size, INT *off, INT *handle)
{
  void *p;
  long offset, realsize;

_IF(cray)
  /* for true cray systems word length computes to 1, hardwire to 8 */
  realsize = 8;
_ELSE
  realsize =  (long) q2 - (long) q1;
_ENDIF

  p = (void *) malloc( *size * realsize );

  if(!p){
    *off = 0;
    return;
  }

  offset = (long ) p - (long) q1;

  /*printf("malloc: allocate %ld \n", p); */

  *handle = (INT) ( (long) p - (long) q1 );

_IF(cray)
 /* cray offset computes directly in words */
  realsize=1;
_ENDIF

  *off = 1 + offset / realsize;  /* add 1 so that Q(*ioff) means what we want */

  /*  printf("check allignment.. \n p=%ld \n q1=%ld \n offset=%ld \n *off = %ld \n q1+*off*8= %ld\n",p,q1,offset,*off, (long) q1 + 8 * *off);*/

}
_IF(64bitpointers)
void freec_(void *q1, long *handle, INT *stat)
_ELSE
void freec_(void *q1, INT *handle, INT *stat)
_ENDIF
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

_ENDIF


_IFN(cray,t3d)
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

_IF(GIGA_DUMP)
void pack3cc_(unsigned long long *a,
              INT *i,
              INT *j,
              INT *k)
{
  *a  = (unsigned long long) *i;
  *a  = *a << 20;
  *a  = *a | (unsigned long long) *j;
  *a  = *a << 20;
  *a  = *a | (unsigned long long) *k;
}

void upack3cc_(unsigned long long *a,
               INT *i,
               INT *j,
               INT *k)
{
  *i =  (INT) (*a >> 40);
  *j =  (INT) ( (*a >> 20)  & 0xFFFFF) ;
  *k =  (INT) (*a & 0xFFFFF);
}
_ELSE
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
_ENDIF
_ENDIF


_IF(rs6000)
_IF(tcgmsg)
/* 
   Adapted error handler. Calls here replace the standard TCGMSG 
   routines in certain cases. Purpose is to clarify meaning of the 
   code and avoid generating a system error message when the error was
   generated internally to GAMESS-UK.
   */

#ifdef CRAYXX
#define pexitc_ PEXITC
#endif

#ifdef NOC_
#define pexitc_ pexitc
#endif

void pexitc_(code)
     INT *code;
{
  (void) fflush(stdout);
  (void) fflush(stderr);
  /*mpc_stopall(code);*/
  abort();
  /*exit(code);*/
}
_ENDIF
_ENDIF

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
_IF(linux,hitachi)

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
_ENDIF
_IF(linux)
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
_ENDIF


/* Added for alex turners matrix code */
void flushout_() {
  fflush(stdout);
}

_IF(signals)

#ifdef CRAYXX
#define iflsh6_ IFLSH6
#define idumpt_ IDUMPT
#define catchcc_ CATCHCC
_IF(sv1)
#define icpu_ ICPU
_ENDIF
#endif

#ifdef NOC_
#define iflsh6_ iflsh6
#define idumpt_ idumpt
#define catchcc_ catchcc
#endif

#include <stdio.h>
#include <sys/types.h>
_IF(sgi)
#include <sys/signal.h>
#include <sys/siginfo.h>
_ELSE
#include <signal.h>
_ENDIF

/*
 *  Handler for USR1 and USR2 signals.
 */

void
_IF(sv1)
sigusr_handler(int sig)
_ELSE
sigusr_handler(int sig, siginfo_t *sip, void *extra)
_ENDIF
{
     switch(sig) {
     case SIGUSR1:
       iflsh6_();
       break;
     case SIGUSR2:
       idumpt_();
       break;
     case SIGXCPU:
       icpu_();
       break;
     }
}

/*
 *  Set up handlers for SIGUSR1 and SIGUSR2
 */

_IF(sv1)
void catchcc_()
  {
    signal(SIGUSR1,sigusr_handler);
    signal(SIGUSR2,sigusr_handler);
    signal(SIGXCPU,sigusr_handler);
  }
_ELSE

void
catchcc_()
{
   struct sigaction sa;
   int ret;
   sa.sa_sigaction = sigusr_handler;
   sigemptyset(&sa.sa_mask);
   sa.sa_flags = SA_SIGINFO;
   ret = sigaction(SIGUSR1, &sa, 0);
   ret = ret || sigaction(SIGUSR2, &sa, 0);
_IF(sgi)
   ret = ret || sigaction(SIGXCPU, &sa, 0);
_ENDIF
   if (ret) {
      perror("sigaction");
   }
}
_ENDIF
_ENDIF

_IF(taskfarm)

#ifdef NOC_
#define rungamessc_  rungamessc
#define taskerrc_   taskerrc
#define setenvc_   setenvc
#define gamess_   gamess
#endif

#ifdef CRAYXX
#define rungamessc_  RUNGAMESSC
#define taskerrc_   TASKERRC
#define setenvc_   SETENVC
#define gamess_   GAMESS
#endif

/* ===========  Function Prototypes/declarations ====================*/
void rungamessc_(INT *init, INT *gamret);
void taskerrc_(INT *gamret);
void setenvc_(STRING_ARGS1,STRING_ARGS2);
/* ===========  end Function Prototypes ====================*/

jmp_buf rungamess_buffer;

void rungamessc_(INT *init, INT *gamret)
{
  int iret;
  /* pass control to fortran code this call will return 0 
        while the buffer is being set up
     it will return again if the jump is executed, this time
        with the code passed to longjmp
  */
  iret = setjmp (rungamess_buffer);
  /* printf("in rungamessc with iret %d\n",iret); */
  
  if(iret){
    *gamret = (INT) iret;
  }else{
    /* printf("in rungamessc with iret %d\n",iret); */
    gamess_(init,&iret);
    *gamret = iret;
  }
}

void taskerrc_(INT *gamret)
{
  INT gerr;
  /* printf("taskerrc recieved *gamret %d\n",*gamret); */
  gerr = (INT) *gamret;
  longjmp(rungamess_buffer,gerr);
}


/* Set the environment variables (used to name ed's 2,3 & 7) */
void setenvc_(STRING_ARGS1,STRING_ARGS2)
{
  char *s, *v;
  int ovwrt,rtn;
  ovwrt=1;
  s = GET_STRING1;
  v = GET_STRING2;
  rtn=setenv(s,v,ovwrt);
}

_ENDIF
_IF(syslog)
#include <syslog.h>
/* Write kernel log message to track errors */
void logsys_(char *txt) {
	int me;
	openlog("GAMESS-UK",LOG_NOWAIT|LOG_PID,LOG_USER);
	me=ipg_nodeid_();
	syslog(LOG_NOTICE,"cpu %d : %s\n",me,txt);
	closelog();
}
_ENDIF

_IF(glibc)
#include <signal.h>
_ENDIF
void segfaulthandler(int sig) {
void *trace[100];
_IF(glibc)
printf("GAMESS segfaulthandler calls backtrace\nrun gdb and x address\n");
backtrace_symbols_fd(trace,backtrace(trace,100),2);
_ENDIF
exit(1);
}

void siginthandler(int sig) {
void *trace[100];
_IF(glibc)
printf("GAMESS siginthandler calls backtrace\nrun gdb and x address\n");
backtrace_symbols_fd(trace,backtrace(trace,100),2);
_ENDIF
exit(1);
}

/* After this segfaulthandler will be called */
void installsegfaulthandler_() {
_IF(glibc)
struct sigaction act;
_IF(linux)
act.sa_handler=(sig_t)segfaulthandler;
_ELSE
act.sa_handler=(void *)&segfaulthandler;
_ENDIF
sigemptyset(&(act.sa_mask));
sigaction(SIGSEGV,&act,NULL);
_ENDIF
}

/* After this the default segv handler will be called */
void uninstallsegfaulthandler_() {
_IF(glibc)
struct sigaction act;
act.sa_handler = SIG_DFL;
sigemptyset (&act.sa_mask);
act.sa_flags = 0;
sigaction(SIGSEGV,&act,NULL);
_ENDIF
}

/* Create a stack trace by Interrupting (ctrl-c) a hanging job */
void installctrlchandler_() {
_IF(glibc)
struct sigaction act;
act.sa_handler=(void *)&siginthandler;
sigemptyset(&(act.sa_mask));
sigaction(SIGINT,&act,NULL);
_ENDIF
}

/* reinstate default ctrl-c */
void uninstallctrlchandler_() {
_IF(glibc)
struct sigaction act;
act.sa_handler = SIG_DFL;
sigemptyset (&act.sa_mask);
act.sa_flags = 0;
sigaction(SIGINT,&act,NULL);
_ENDIF
}

void fgetdate_(char zdate[])
{
    time_t timet;

    timet=time(NULL);
    strncpy(zdate,ctime(&timet),24);

}
