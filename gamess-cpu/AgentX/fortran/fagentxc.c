
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

#include<agentx.h>
#include<fagentxInternals.h>

#include<stdlib.h>
#include<string.h>

void ax_fortran_string(char* result, unsigned int resultlen,char* retval){

  unsigned int i=0;

  if(!retval)goto cleanup;
  if(!result)goto cleanup;
  if(!resultlen)goto cleanup;

  for(i=0;i<resultlen;i++){
    if(i<strlen(retval)){
      *(result+i)=*(retval+i);
    }
    else{
      *(result+i)=' ';
    }
  }
 cleanup:
  return;
}

char *ax_null_terminate(char *s, int n)
{
  /* remove trailing white space and null terminate
     a Fortran string */
  int i;
  char *t;
  for(i = n - 1; *(s+i) == ' ' && i; i--)
    ;
  t = (char *) malloc(sizeof(char) * (i+2));
  if(!t){t=NULL;goto cleanup;}
  strncpy(t,s,i+1);
  *(t+i+1) = '\0';
 cleanup:
  return t;
}

fortint axcwrapproxy(fortintc* proxylen, char* proxy){
  int retval=0;
  char *t = ax_null_terminate(proxy,(int)(*proxylen));
  retval=axProxy(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapgeturi(fortintc* urilen, char* uri){
  int retval=0;
  char *t = ax_null_terminate(uri,(int)(*urilen));
  retval=axGetUri(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapdatageturi(fortintc* urilen, char* uri){
  int retval=0;
  char *t = ax_null_terminate(uri,(int)(*urilen));
  retval=axDataGetUri(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapselect(fortintc* entitylen, char* entity){
  int retval=0;
  char *t = ax_null_terminate(entity,(int)(*entitylen));
  retval=axSelect(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapselectnext(){
  int retval=0;
  retval=axSelectNext();
  return (fortint)retval;
}

fortint axcwrapdeselect(){
  int retval=0;
  retval=axDeselect();
  return (fortint)retval;
}

fortint axcwrapbaseuri(fortintc* local_base_uri_stringlen, char* local_base_uri_string){
  int retval=0;
  char *t = ax_null_terminate(local_base_uri_string,(int)(*local_base_uri_stringlen));
  retval=axBaseUri(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapparserfinish(){
  int retval=0;
  retval=axParserFinish();
  return(fortint)retval;
}

fortint axcwrapparserstart(){
  int retval=0;
  retval=axParserStart();
  return(fortint)retval;
}

fortint axcwraprefine(fortintc* expressionlen, const char* expression){
  int retval=0;
  char *t = ax_null_terminate((char*)expression,(int)(*expressionlen));
  retval=axRefine(t);
  free(t);
  return (fortint)retval;
}

fortint axcwrapformat(fortint *dataFormat){
  int retval=0;
  retval=axFormat((int)(*dataFormat));
  return(fortint)retval;
}
  
fortint axcwrapcurrent(){
  int retval=0;
  retval=axCurrent();
  return (fortint)retval;
}

fortint axcwraplen(){
  int retval=0;
  retval=axLen();
  return (fortint)retval;
}

void axcwrapvalue(fortintc* strvarlen, char* strvar){
  char* retval=NULL;
  retval=axValue();
  ax_fortran_string(strvar,(unsigned int)(*strvarlen),retval);
}

fortint axcwrapbuffer(fortint* size){
  int retval=0;
  retval=axBuffer((int)(*size));
  return(fortint)retval;
}

fortint axcwrapnamespace(fortintc* namespacefURIlen, char* namespacefURI){
  int retval=0;
  char *t = ax_null_terminate(namespacefURI,(int)(*namespacefURIlen));
  retval=axNamespace(t);
  free(t);
  return(fortint)retval;
}

fortint axcwrapcache(fortint *cachesize){
  int retval=0;
  axCache((int)(*cachesize));
  return(fortint)retval;
}

fortint axcwrapselectno(fortint *no){
  int retval=0;
  axSelectNo(*no);
  return(fortint)retval;
}
