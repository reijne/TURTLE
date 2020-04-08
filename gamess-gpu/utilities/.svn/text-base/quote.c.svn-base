#include <stdio.h>
#include <string.h>

/* change " into ' if " is not between quotes; ('"') is not changed */

#define MAXCHARS 90

main()
{
   char regel[MAXCHARS];
   int i;
   int ilen,ic;
   const int mask=1;

   while(fgets(regel,MAXCHARS,stdin)!=NULL){
     ilen=strlen(regel);
     if ((regel[0]!='c')&&(regel[0]!='*')&&(regel[0]!='C')&&(ilen>6)){
        ilen=strlen(regel);
        if(regel[5]==' ')ic=0;
        for(i=0;i<ilen;i++){
           if (regel[i]=='\'') ic^=mask;
           if ((regel[i]=='\"')&&!ic) regel[i]='\'';
        }
     }
     fputs(regel,stdout);
     fflush(NULL);
   }
}            
               


