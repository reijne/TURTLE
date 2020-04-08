#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define New(t)          (t *)calloc(1, sizeof(t))
#define Chkp(p) if(!(p)){fprintf(stderr,"Memory allocation error\n"); exit(1);}
/*
   simple code converter

   This version is adapted (July 94) to allow insertion of code
   fragments

   Adapted (April 1995 to allow "fake" directory specification
   via CCNVRTDIR.

   Jan 02 - added ability to add an error trap in the source

*/

void windup(char *prog);
void output(char *prog, char *s);
void strcpym(char **s1, char *s2);
/* renamed getline to psh_getline to avoid name clash on Windows */
int psh_getline(FILE *fp, char *line, int max);
int matchtag(char *s);
int split_line(char *line, char *strings, int maxs, int maxc);

#define LENBUFF 4096
char buf1[LENBUFF], buf2[LENBUFF + 10];
char leader[10], key[LENBUFF], trailer[10];

char filename[4096],leaddir[4096]; /* directory and filename */
char *tmpdir;

char *file, *outfile;
int fflag = 0;
int  nl,nt,i1,i2;
int debug=0;
int no;
int ibl;
FILE *fpin, *fpout, *fopen();

struct flag_struct { char *tag; struct flag_struct *prev;} *new, *last, *flag;

#define PEND 0
#define COPY 1
#define SKIPINCLUDE -3
#define OVER -2

#define IF 3
#define ELSE 2
#define ELSEIF 4
#define END 5
#define INCLUDE 6
#define ENDINCLUDE 7
#define ERROR 8

struct stack_struct { int status; struct stack_struct *prev;} *top, *newtop;

main(int argc,
     char *argv[])
{
  int i, islash,h1,h2,iret,shift,inside,copy,type,clean,control;
  int com();
  char *strst();
  char *getenv();
  char *c;
  FILE *fp[20];
  int ilevel;

  last = NULL;
  top = NULL;
  inside = 0;
  clean = 0;
  no = 0;
  ilevel = 0;

  /* parse args */

  for(i=1;i<argc;i++){
    if(!strcmp(argv[i],"-f")){        /* specify input/output file */
      strcpym(&file,argv[++i]);
      fflag = 1;
    }else if(!strcmp(argv[i],"-c")){  /* produce clean code */
      clean = 1;
    }else if(!strcmp(argv[i],"-d")){  /* produce clean code */
      debug = 1;
    }else if(!strcmp(argv[i],"-n")){  /* no output - just record changes */
      no = 1;
    }else if(argv[i][0] == '-'){ 
      fprintf(stderr,"%s: bad flag %s\n",argv[0],argv[i]);
      exit(1);
    }else{          /* assume it is a code flag */
      new = New(struct flag_struct);
      new->prev = last;
      last = new;
      strcpym(&(new->tag),argv[i]);
    }
  }

  /* CCNVRTDIR - leading directory specification, used for
     include files */

  strcpy(leaddir,".");
  tmpdir = getenv("CCNVRTDIR");
  if(tmpdir)strcpy(leaddir,tmpdir);


  /* dont let the the user trash the control data */

  if(clean && fflag){
    fprintf(stderr,"%s: -c and -f flags are not compatible\n",argv[0]);
    exit(1);
  }

/* set up io streams */

  if(fflag){
    fpin = fopen(file,"r");
    if(!fpin){
      fprintf(stderr,"%s: attempt to open input file %s failed\n",argv[0],file);
      exit(1);
    }
    for(islash = -1,i=0; i < strlen(file); i++)if(file[i] == '/')islash=i;
    strcpy(buf1,"/tmp/");
    if(islash != -1)islash++;
    else islash=0;
    strcat(buf1,file+islash);
    strcpym(&outfile,buf1);
    fpout = fopen(outfile,"w");
    if(!fpout){
      fprintf(stderr,"%s: attempt to open output file %s failed\n",argv[0],outfile);
      exit(1);
    }
  }else{
    fpin = stdin;
    fpout = stdout;
  }

  if(debug){
    printf("flags set :");
    for(flag=last;flag;flag=flag->prev)printf("%s ",flag->tag);
    printf("\n");
  }

  /* process file */

  copy = 1;
  iret = 0;
  while(1){
  xxxx:
    iret = psh_getline(fpin,buf1,LENBUFF);
    if(iret == -1){
      fprintf(stderr,"%s: line too long\n",argv[0]);
      exit(1);
    }

  /* done ? */

    if(iret == -2){
      if(ilevel == 0){
	windup(argv[0]);
      }else{
	if(debug)printf("include finished at level %d\n",ilevel);
	/* we have finished parseing an included file */
	ilevel--;
	/* restore pointer */
	fpin = fp[ilevel];

	/* change state to SKIPINCLUDE */
	newtop = New (struct stack_struct);
	newtop->prev = top;
	top = newtop;
	top->status = SKIPINCLUDE;
	inside = 1;

	/* try and load next record of old file */
	goto xxxx;
      }
    }
 /* 
   check whether it is a control statment

   [comment leader]#include file#[comment trailer]  
   [comment leader]#endinclude file#[comment trailer]  
   [comment leader]#if flaglist#[comment trailer]  
   [comment leader]#elseif flaglist#[comment trailer]  
   [comment leader]#else#[comment trailer]  
   [comment leader]#endif#[comment trailer]    (#end# also allowed)
   [comment leader]#error#[comment trailer]

    (no leader -> leader = #)
 */

    control = 0;
    h1 = -1;
    h2 = 0;

    for(c=buf1;*c;c++){
      if(*c == '#'){
	if(h1 == -1)h1 = c - buf1;
	else if (!h2)h2 = c - buf1;
	else {
	  h1 = h2;
	  h2 = c - buf1;
	}
      }
    }
    if(h2){    /* try and parse */
      for(i=0; i < h2 - h1 -1; i++)key[i] = buf1[h1 + i + 1];
      key[i] = '\0';

      control=1;
      if(strst(key,"endinclude ") == key){
	type = ENDINCLUDE; 
	shift = 11;
      }else if(strst(key,"end") == key){
	type = END;
	shift = 4;
      }else if(strst(key,"endif") == key){
	type = END;
	shift = 6;
      }else if(strst(key,"elseif ") == key){
	type = ELSEIF;
	shift = 7;
      }else if(strst(key,"else") == key){
	type = ELSE;
	shift = 5;
      }else if(strst(key,"if ") == key){
	type = IF; 
	shift = 3;
      }else if(strst(key,"include ") == key){
	type = INCLUDE; 
	shift = 8;
      }else if(strst(key,"error ") == key){
	type = ERROR; 
	shift = 6;
      }else{
	control=0;
      }
    }

    if(control){
      if ( h1 > 9 ){
	fprintf(stderr,"%s: comment header is too long on line %s\n",argv[0],buf1);
	exit(1);
      }
      for(i = 0; i < h1; i++)leader[i] = buf1[i];
      leader[i] = '\0';
      nl = i;
      if(nl == 0){  /* default comment leader */
	nl = 1;
	strcpy(leader,"#");
      }

/* copy trailer (if it is nonblank) */

      ibl = 1;
      for(i = 0; buf1[h2 + i + 1]; i++)if(buf1[h2 + i + 1] != ' ')ibl = 0;

      if(ibl){
	nt = 0;
      }else{
	for(i = 0; buf1[h2 + i + 1]; i++){
	  if ( i > 9 ){
	    fprintf(stderr,"%s: comment trailer is too long %s\n",argv[0],buf1+h2+1);
	    exit(1);
	  }
	  trailer[i] = buf1[h2 + i + 1];
	}
	trailer[i] = '\0';
	nt = i;
      }

      if(debug)printf("leader = %s, key = %s, trailer = %s\n",leader,key,trailer);

      if(shift){
	for(i=shift;key[i] == ' ';i++);
	strcpy(buf2,key+i);
	strcpy(key,buf2);

      }
/* 
    apply control logic 
*/

      /* suppress effect of control statments for SKIPINCLUDE region */
      if(top && top->status == SKIPINCLUDE){
	switch(type){
	case INCLUDE:
	case ENDINCLUDE:
	  break;
	default:
	  type =  0;
	  break;
	}
      } 

      switch (type){
      case IF:  /* push stack */
	newtop = New (struct stack_struct);
	newtop->prev = top;
	top = newtop;
	if(inside && !copy){
	  top->status = OVER;
	}else{
	  top->status = matchtag(key);
	}
	break;
      case ELSEIF:
	if(!top){
	  fprintf(stderr,"%s: elseif without if %s\n",argv[0],buf1);
	  exit(1);
	}
	switch(top->status){
	case COPY:
	  top->status = OVER;
	  break;
	case PEND:
	  top->status = matchtag(key);
	  break;
	case OVER: /* null case */
	  break;
	}
	break;
      case ELSE:
	if(!top){
	  fprintf(stderr,"%s: else without if %s\n",argv[0],buf1);
	  exit(1);
	}
	switch(top->status){
	case COPY:
	  top->status = OVER;
	  break;
	case PEND:
	  top->status = COPY;
	  break;
	case OVER: /* null case */
	  break;
	}
	break;
      case END:  /* pop stack */
	if(!top){
	  fprintf(stderr,"%s: unmatched %s\n",argv[0],buf1);
	  exit(1);
	}
	newtop=top->prev;
	free(top);
	top = newtop;
	/* restore status */
	break;
      case INCLUDE:

	if(top && top->status == SKIPINCLUDE){
	  /* push another level on so we can count the endincludes */
	  newtop = New (struct stack_struct);
	  newtop->prev = top;
	  top = newtop;
	  top->status = SKIPINCLUDE;
	} else {
	  /* save old file pointer */
	  fp[ilevel]=fpin;
	  /* create new file pointer */ 
	  ilevel++;
	  if(ilevel >= 20){
	    fprintf(stderr,"too many include levels (>19) - check for recursive includes\n");
	    exit(1);
	  }
	  /* note that the stack is unchanged by the inclusion,
	     until the eof of the included file is reached 
	     when SKIPINCLUDE is used */

	  strcpy(filename,leaddir);
	  strcat(filename,"/");
	  strcat(filename,key);

	  fpin = fopen(filename,"r");
	  if(!fpin){
	    fprintf(stderr,"failed to include file %s\n",filename);
	    exit(1);
	  }
	}
	break;
      case ENDINCLUDE:
	if(!top){
	  fprintf(stderr,"%s: unmatched %s\n",argv[0],buf1);
	  exit(1);
	}
	switch(top->status){
	case SKIPINCLUDE:
	  newtop=top->prev;
	  free(top);
	  top = newtop;
	}
	break;

      case ERROR:  /* push stack */

	if(!top || top->status == COPY){
	  fprintf(stdout,"ccnvrt: %s\n",key);
	  fprintf(stderr,"ccnvrt: %s\n",key);
	  exit(-1);
	}
	break;

      case 0:
	break;
      default:
	fprintf(stderr,"unknown control %s\n",buf1);
	exit(1);
      }
      if(top){
	copy = (top->status > 0);
	inside = 1;
      }else{
	inside = 0;
      }

      /* output the control record */
      if(debug && ! clean)printf("control:");
      if(!clean && (!top || (top->status != SKIPINCLUDE)))output(argv[0],buf1);
      else if(no)no++;

    }else{
      if(!inside){
	/* no switches are active, save time by copying directly */
	if(debug)printf("plain:");
	output(argv[0],buf1);
      }else if (top->status == SKIPINCLUDE){
	/* 
	   dont even comment ... as this piece of the
	   the file has been replaced by an inclusion
	*/
      }else if (!copy){
	if(!clean){ 
	  if(com()){
	    if(debug)printf("copying comment:");
	    output(argv[0],buf1);
	  }else{
	    if(no)no++;
	    if(debug)printf("commenting:");
	    sprintf(buf2,"%s%s%s",leader,buf1,trailer);
	    output(argv[0],buf2);
	  }
	}else{
	  if(no)no++;
	}
      } else{
	if(com()){
	  if(no)no++;
	  if(debug)printf("uncommenting:");
	  if(i2)buf1[i2] = '\0';
	  output(argv[0],buf1+i1);
	}else{
	  if(debug)printf("copying noncomment:");
	  output(argv[0],buf1);
	}
      }
    }
  }
}
/* see if the current line is already a comment */
int com()
{
  int i;
  char *strst();
  if(nl && nt){   /* both header and leader */
    i1 = strst(buf1,leader) - buf1;
    if(i1 < 0)return 0;
    for(i=0;i<i1;i++)if(buf1[i] != ' ')return 0;
    i2 = strst(buf1 + i1 + 1,trailer) - buf1;
    if(i2 < 0)return 0;
    for(i=i2+nt; i<strlen(buf1); i++)if(buf1[i] != ' ')return 0;
    i1 = i1 + nl;
    i2 = i2 - 1;
    return 1;
  }else{         /* the string must start with the leader */
    i1 = nl;
    i2 = 0;
    return (strst(buf1,leader) == buf1);
  }
}
/*
    see if any component of s matches any of the
    command line keys
    !sgi and ! sgi are hits if sgi was not set

*/
int matchtag(char *s)
{
  int ns, i, miss, neg;
  char split[20][100];
  ns = split_line(s,split[0],20,100);
  for(i=0;i<ns;i++){
    neg = 0;
    if(!strcmp(split[i],"!")){
      neg = 1;
      i++;
    }else if(split[i][0] == '!'){
      neg = 1;
      strcpy(buf2,split[i]+1);
      strcpy(split[i],buf2);
    }
    if(neg){
      miss = 0;
      for(flag=last;flag;flag=flag->prev)if(!strcmp(flag->tag,split[i]))miss = 1;
      if(!miss)return COPY;
    }else{
      for(flag=last;flag;flag=flag->prev)if(!strcmp(flag->tag,split[i]))return COPY;
    }
  }
  return PEND;
}
/* 
    termination routine
*/
void windup(char *prog)
{
  int iret;
/* check the stack */
  if(top){
    fprintf(stderr,"%s: if block not closed at end of file\n",prog);
    exit(1);
  }
  if(fflag && ! no){  /* try and move back the file */
    fclose(fpout);
    sprintf(buf1,"cp %s %s",outfile,file);
    iret = system(buf1);

    if(iret){
      fprintf(stderr,"%s: attempt to replace file %s failed\n",prog,file);
      exit(1);
    }
    unlink(outfile);
  }
  if(no){
    no--;
    printf("changes = %d\n",no);
  }
  exit(0);
}
/*
    fprintf with an error trap and bypass if -n
*/
void output(char *prog,
	    char *s)
{
  int iret;
  if(no)return;
  iret = fprintf(fpout,"%s\n",s);
  if(iret < 0){
    fprintf(stderr,"%s: output error on %s\n",prog,outfile);
    exit(1);
  }
}
/********************************************************************
 * strst1 : index a single character q in a string p and return posn
 ********************************************************************/
int strst1(char *p, char *q)
{
  char *p1;
  for(p1 = p  ;*p1 != '\0' ;p1++)
               if(*p1 == *q)return p1-p+1;
   return 0;
}
/********************************************************************
 * strst : index the occurence of q within p returning string p
 ********************************************************************/
char *strst(char *p, char *q)
{
  char *p1,*q1,*p2;
  for(p1 = p  ;*p1 != '\0' ;p1++){
    for(p2 = p1, q1 = q; (*p2 == *q1) && (*q1 != '\0') ;q1++,p2++)
      ;
    if(*q1 == '\0')
      return p1;
  }
  return (char *) NULL;
}
void strcpym(char **s1, char *s2)
{
  *s1 = malloc(sizeof(char) * (strlen(s2) + 1));
  strcpy(*s1,s2);
}
/********************************************************************
 * psh_getline()
 * read a line (up to max chars) from a stream fp into line
 * remove trailing \n if present
 * return code - 0 if OK
 * -1 if buffer overflow 
 * -2 if EOF occurs
 ********************************************************************/
int psh_getline(FILE *fp,
	    char *line,
	    int max)
{
  int i;
  char *ret, *fgets();
  
  *line = '\0';
  ret = fgets(line,max,fp);
  if(!ret){
    return -2;
  }

  i = strst1(line,"\n");
  if (i != 0 && i < max)   /* newline lies within string */
    {
      *(line+i-1) = '\0';
      return i-1;    
    }
  else                    /* no newline */
    {
      *(line + max - 1) = '\0';
      return -1;
    }
}
/**************************************************************
 * split a line into strings
 * return the number of strings
 **************************************************************/
int split_line(char *line,
	       char *strings,
	       int maxs,
	       int maxc)
{
  int ns, nc, instring;
  char *p;
  /* parse into fields, = is a field on its own*/
  ns = -1;
  nc = 0;
  instring = 0;
  for(p = line; *p != '\0' && *p != '\n' ;p++)
    {
      if(*p == ' ')
	{
	  if(!instring)
	    ;
	  else
	    {
	      *(strings+maxc*ns+nc) = '\0';
	      nc = 0;
	      instring = 0;
	    }
	}
      else if(*p == '"')
	if(instring)
	  {
	    printf("error \" occurred inside a quoted string");
	    return -1;
	  }
	else
	  {
	    nc = (int) (strst((p+1),"\"")-p-1);
	    strncpy(strings+(++ns)*maxc,p+1,nc);
	    *(strings+ns*maxc+nc) = '\0';
	    p = p+nc+1;
	    nc = 0;
	  }
      else if(*p == '=')
	{
	  if(instring)
	    {
	      *(strings+ns*maxc+nc) = '\0';
	      instring = 0;
	    }
	  *(strings+(++ns)*maxc) = '=';
	  *(strings+ns*maxc) = '\0';
	  nc = 0;
	}   
      else
	{
	  /* other character */
	  if(!instring)
	    {
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
