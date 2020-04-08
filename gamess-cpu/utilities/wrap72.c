#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#define MAXBUFF 132
char buff[MAXBUFF];
char buff72[73];

main(argc, argv)
     int argc;
     char *argv[];
{
  int l;
  char *p;
  int iret;
  strcpy(buff72+72,"");
  while(1){
    iret = psh_getline(stdin,buff,MAXBUFF-1);
    if(iret >= 0){
      if(buff[0] == 'c' || buff[0] == 'C' || buff[0] == '*'){
	printf("%s\n",buff);
      }else{
	l = strlen(buff);
	if(l <= 72){
	  printf("%s\n",buff);
	}else{
	  fprintf(stderr,"wrapped  line: %s\n",buff);
	  strncpy(buff72,buff,72);
	  printf("%s\n",buff72);
	  p = buff + 72;
	  strcpy(buff72,"     +");
	  while(p < buff + l){
	    strncat(buff72,p,66);
	    printf("%s\n",buff72);
	    p += 66;
	  }
	}
      }
    }else if (iret == -1){
      strcpy(buff + MAXBUFF-1,"");
      fprintf(stderr,"%s line too long for buffer: %s\n",argv[0],buff);
      exit(1);
      break;
    }else if (iret == -2){
      exit(0);
    }
  }
}
/********************************************************************
 * psh_getline()
 * read a line (up to max chars) from a stream fp into line
 * remove trailing \n if present
 * return code - 0 if OK
 * -1 if buffer overflow 
 * -2 if EOF occurs
 ********************************************************************/
int psh_getline(fp, line, max)
FILE *fp;
char *line;
int max;
{
  int i;
  char *ret;

  *line = '\0';
  if(!fgets(line,max,fp))
    return -2;

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
