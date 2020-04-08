/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/strdup.c,v 1.1.1.5 2007-10-30 10:14:16 jmht Exp $ */

#include <stdlib.h>
extern char *strcpy();
extern size_t strlen();

char *strdup(s)
    char *s;
{
  char *new;

  if ((new = malloc((size_t) (strlen(s)+1))))
     (void) strcpy(new,s);

  return new;
}
