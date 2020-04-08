#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#define MAXCHK 20
#define MAXC 300
#define MAXS 30

int debug = 1, silent=0;

main(argc,argv)
     int argc;
     char *argv[];
{
  char key[MAXC];
  char buff[MAXC];
  char vtab[MAXC];
  char grepper[MAXCHK][MAXC];
  char label[MAXCHK][MAXC];
  char bits[MAXS][MAXC];
  char refmach[MAXS][MAXC];
  char *strst();
  char *rem_spaces();
  char err_line[MAXC];
  double answer[MAXCHK], tester, test;
  double toler[MAXCHK];
  int time[MAXCHK];
  int match, finished,ngrep[MAXCHK],ncheck[MAXCHK],nskip[MAXCHK],nhits,nchk,ichk,ii,ns;
  int code, retval, err_det, err_cnt;

  FILE *fp, *fopen();
  if(argc != 3){
    fprintf(stderr,"usage: %s jobname logfile\n",argv[0]);
    exit(-1);
  }
/* 
    open validation table
*/
  if(getenv("GAMESS_VTAB"))strcpy(vtab,getenv("GAMESS_VTAB"));
  else strcpy(vtab,"VTAB");
  fp = fopen(vtab,"r");
  if(!fp){
    fprintf(stderr,"cant open %s\n",vtab);
    printf("%s cant verify - no table\n",argv[1]);
    exit(-1);
  }

  if(getenv("GAMESS_VALIDATE_QUIETLY")){
    silent=1;
  }

/*
    obtain check information
*/
  finished = match = nchk = 0;


  while (! finished){
    if(fscanf(fp,"%s",key) == 1){
      psh_getline(fp,buff,MAXC);
      if(!strcmp(key,argv[1])){
	time[nchk]=0;
	ns = split_line(buff,bits,MAXS,MAXC);
	if( ns != 7 && ns != 8){
	  printf("%s  test %d bad - wrong number of fields in table \n",argv[1],nchk+1);
	  exit(-1);
	}
	strcpy(label[nchk],bits[0]);
	strcpy(grepper[nchk],bits[1]);
        nskip[nchk] = atoi(bits[2]);
	ngrep[nchk] = atoi(bits[3]);
	ncheck[nchk] = atoi(bits[4]);
	if(sscanf(bits[5],"%lf",&answer[nchk]) != 1){
	  printf("%s  test %d bad - target is not a number\n",argv[1],nchk+1);
	  exit(-1);
	}
	if(!strcmp(bits[6],"time")){
	  time[nchk]=1;
	  if(ns == 8)sprintf(refmach[nchk],"* %s",bits[7]);
	  else strcpy(refmach[nchk],"");
	}else if(sscanf(bits[6],"%lf",&toler[nchk]) != 1){
	  printf("%s  test %d bad - tolerance percentage is not a number or the string \"time\" \n",argv[1],nchk+1);
	  exit(-1);
	}
	nchk++;
	match = 1;
      }
    }else{
      finished = 1;
    }
  }
  if(! match){
    fprintf(stderr," no validation table entry for %s\n",argv[1]);
    printf("%8s cant verify - not in table\n",argv[1]);
    exit(-2);
  }

/* 
    check the file
*/
  fp = fopen(argv[2],"r");
  if(!fp){
    fprintf(stderr,"cant open %s\n",argv[2]);
    printf("%8s cant verify - cant find file %s\n",argv[1],argv[2]);
    exit(-1);
  }
  retval= 0;
  err_det=0;
  finished = match = nhits = ichk = 0;
  while (! finished){
    if(psh_getline(fp,buff,MAXC) >= 0){
      if(strst(buff,"error detected")) {
        err_det=1;
	err_cnt=0;
      }
      if(err_det) {
        err_cnt++;
	if (err_cnt==4) {
	  strcpy(err_line,rem_spaces(buff));
          if (strst(err_line,"i/o error; logical file =ed0")) {
            strcpy(err_line,"ecp library not provided (ed0)");
          }
          if (strst(err_line,"i/o error; file = table")) {
            strcpy(err_line,"file TABLE not provided");
          }
	  printf("%20s (%30s) %s\n",argv[1],label[ichk],err_line);
	  err_det=0;
	}
	retval = 4;
      }      
      if(ichk < nchk && strst(buff,grepper[ichk])){
	nhits++;
	if(nhits == ngrep[ichk]){
	  if(ncheck[ichk] == 0){
	    /* Simple existence check */
	    if(!silent)printf("%20s (%30s) test %d OK \n",argv[1],label[ichk],ichk+1);
	  }else{
	    /* Numerical check */
	    for(ii=0;ii<nskip[ichk];ii++)psh_getline(fp,buff,MAXC);
	    ns = split_line(buff,bits,MAXS,MAXC);
	    if(sscanf(bits[ncheck[ichk]-1],"%lf",&tester) != 1){
	      printf("%20s (%30s) test %d failed - target was not a number\n",argv[1],label[ichk],ichk+1);
	      if(debug)for(ii=0;ii<ns;ii++)fprintf(stderr," %s \n",bits[ii]);
	      exit(-1);
	    }
	    if(time[ichk]){
	      /* The field is a reference time, don't bother checking for validity */ 
	      test = answer[ichk] / tester;
	      printf("%20s (%30s) test %d Speed %f %s\n",argv[1],label[ichk],ichk+1,test,refmach[ichk]);
	    }else{
	      test = 100 * fabs( (double) 1.0 - tester/answer[ichk]);
	      if(test < toler[ichk]){
		if(!silent)printf("%20s (%30s) test %d OK \n",argv[1],label[ichk],ichk+1);
	      }else{
		printf("%20s (%30s) test %d error %lf percent  %lf %lf \n",argv[1],label[ichk],ichk+1,test,tester,answer[ichk]);
		retval = 2;
	      }
	    }
	  }
	  nhits = 0;
	  ichk ++;
	}
	if(ichk == nchk+1)finished=1;
      }
    }else{
      finished = 1;
    }
  }
  if(ichk != nchk){
    printf("%20s (%30s) failed - missing result \n",argv[1],label[ichk]);
    exit(-1);
  }

  /* finally check return code if present appended */
  finished=0;
  while (!finished) {
    if(psh_getline(fp,buff,MAXC) >= 0){
      if(strst(buff,"Return Code")){
	ns = split_line(buff,bits,MAXS,MAXC);
	if(sscanf(bits[3],"%d",&code) != 1){
	  printf("failed to locate return code file=%s \n",argv[1]);
	  exit(-1);
	}
	if(code){
	  printf("%8s  *************  Return Code =  %d \n",argv[1],code);
	  exit(3);
	}
      }
    }else{
      finished=1;
    }
  }
  exit(retval);
}
/********************************************************************
 * psh_getline()
 * read a line (up to max chars) from a stream fp into line
 * remove trailing \n if present
 * return code - 0 if OK
 * -1 if EOF occurs
 ********************************************************************/
int psh_getline(fp, line, max)
FILE *fp;
char *line;
int max;
{
  int i=0;

  do {
      *line = '\0';
      if(!fgets(line,max,fp))
        return -1;
      i = strst1(line,">");
     } while (i==2);

  i = strst1(line,"\n");
  if (i != 0 && i < max)   /* newline lies within string */
    {
      *(line+i-1) = '\0';
      return i-1;    
    }
  return 0;
}
/**************************************************************
 * split a line into strings
 * return the number of strings
 **************************************************************/
int split_line(line,strings,maxs, maxc)
char *line;
char *strings;
int maxs,maxc;
{
  int ns, nc, instring;
  char *p;
  char * strst();
  /* parse into fields */
  ns = -1;
  nc = 0;
  instring = 0;
  for(p = line; *p != '\0' && *p != '\n' ;p++){
    if(*p == ' '){
      if(!instring)
	;
      else {
	*(strings+maxc*ns+nc) = '\0';
	nc = 0;
	instring = 0;
      }
    } else if(*p == '"')
      if(instring) {
	printf("error \" occurred inside a quoted string");
	return -1;
      } else {
	nc = (int) (strst((p+1),"\"") - p - 1);
	strncpy(strings+(++ns)*maxc,p+1,nc);
	*(strings+ns*maxc+nc) = '\0';
	p = p+nc+1;
	nc = 0;
      } else{
	/* other character */
	if(!instring) {
	  if(ns > maxs - 2)return -2; /* too many strings to store */
	  instring = 1;
	  ns++;
	}
	if(nc < maxc -1)
	  *(strings+ns*maxc+nc++) = *p;
	else
	  printf("WARNING string truncated to %d characters\n",maxc);
      }
  } /* end of loop over characters */
  if(instring)
    *(strings+ns*maxc+nc) = '\0';
  
  ns++;
/*
    for(i = 0;i <  ns; i++)
      printf("string %s\n",strings + i*maxc);
*/
  return ns;
}
/********************************************************************
 * strst1 : index a single character q in a string p and return posn
 ********************************************************************/
int strst1(p,q)
char *p, *q;
{
  char *p1;
  for(p1 = p  ;*p1 != '\0' ;p1++)
               if(*p1 == *q)return p1-p+1;
   return 0;
}
/********************************************************************
 * strst : index the occurence of q within p returning string p
 ********************************************************************/
char *strst(p,q)
char *p, *q;
{
  char *p1,*q1,*p2;
  for(p1 = p  ;*p1 != '\0' ;p1++)
    {
      for(p2 = p1, q1 = q; (*p2 == *q1) && (*q1 != '\0') ;q1++,p2++)
	;
      if(*q1 == '\0')
	return p1;
    }
  return NULL;
}
/********************************************************************
 * rem_spaces : removes leading and trailing spaces from p and
 *              store them in q                                      
 ********************************************************************/
char *rem_spaces(p)
char *p;
{
  int first;
  char *p1;
  char *q;
  int len, still_space;
  first=1;
  q=p;
  for(p1 = p;*p1!= '\0';p1++) {
    if ((*p1 != ' ') && (first)) {
      q=p1;
      first=0;
    }
  }
  if (q==p) q=p1;
  len=strlen(q);  
  still_space=1;
  for (p1=q+(len-1) ; p1 >= q; p1--) {
    if ((*p1!=' ') && (still_space)) {
      still_space=0;
      *(p1+1)='\0';
    }
  }
  return q;
}
#ifdef NOSTRLEN
/********************************************************************
 * strlen : Returns the length of specified string                   
 ********************************************************************/
int strlen(p)
char *p;
{
  char *p1;
  int ii;
  ii=0;
  for (p1=p;*p1 !='\0' ;p1++) {
    ii++;
  }
  return ii;
}
#endif
