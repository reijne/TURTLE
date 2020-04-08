
#include<agentxInternals.h>
#include<agentx.h>
#include<axArray.h>
#include<dictionary.h>
#include<iniparser.h>

#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

/* static function prototypes */
static int getINIFname( char *, size_t );

/* initialise statics */

static int axStateRec = 0;
static xmlChar **charData = NULL;
static xmlParserCtxtPtr pctxt = NULL;
static int cachehits = 0;
static int useCache = 1;
static int maxEvalCacheSize = 3000;
static int currentEvalCacheSize = 0;
static int linksResolved = 0;
static int charbuffer = -1;
static int dataFormat = 1;
static int supdocs = 0;
static xmlXPathContextPtr xpathCtx = NULL;
static xmlChar* base_uri_string = NULL;
static xmlChar* proxyuri = NULL;

static struct ctxStr* ctxStrPtr = NULL;
static struct ctxStr* ctxStrOrigin = NULL;

static struct cacheBtree* xptrCacheOrigin = NULL;
static struct cacheBtree* propertyCacheOrigin = NULL;
static struct cacheBtree* classCacheOrigin = NULL;
static struct cacheBtree* typeCacheOrigin = NULL;

static xmlChar* namespaceURI = NULL;

int axCheckState(int state, int expState){

  return ( (axStateRec & state) == expState );

}

int axParserStart(){

  /* initialise parser */

  if ( axCheckState ( AX_INIT, AX_ON ) ) {
#ifdef DEBUG
    fprintf(stderr, "Error: parser in wrong state: axStateRec = %d\n", axStateRec);
#endif
    return -1;
  }

  xmlInitParser();

  base_uri_string = (xmlChar*) malloc ( ( MAX_NS_LEN + 1 ) * sizeof (xmlChar) );
  strcpy ( (char*) base_uri_string, "http://www.grids.ac.uk/eccp/owl-ontologies#" );

  namespaceURI = (xmlChar*) malloc ( ( MAX_NS_LEN + 1 ) * sizeof (xmlChar) );
  strcpy ( (char*) namespaceURI, "http://www.grids.ac.uk/eccp/ns#");

  axStateRec += AX_INIT;

  xptrCacheOrigin = NULL;
  propertyCacheOrigin = NULL;
  classCacheOrigin = NULL;
  typeCacheOrigin = NULL;

  return(0);
}

int axFormat(int dataFormatVar){

  dataFormat = dataFormatVar;

  return(0);
}

int axDataGetUri(char* uri){

  if ( dataFormat == 1 ){
    return xml_get_uri( (xmlChar*) uri );
  }
  return -1 ;
}

int axSelect(char* term){

  if ( axCheckState( AX_INIT, AX_OFF ) ) {
#ifdef DEBUG
    fprintf(stderr, "Error: parser in wrong state: axStateRec = %d\n", axStateRec);
#endif
    return -1;
  }

  if( dataFormat == 1 ){
    return xml_select( (xmlChar*)term, 0 );
  }
  return -1;
}

int axSelectType(char* term,int type){

  if ( axCheckState( AX_INIT, AX_OFF ) ) {
#ifdef DEBUG
    fprintf(stderr, "Error: parser in wrong state: axStateRec = %d\n", axStateRec);
#endif
    return -1;
  }

  if ( dataFormat == 1 ){
    return xml_select( (xmlChar*)term, type );
  } 
  return -1;
}

int axSelectNext(){

  if(dataFormat==1){
    return xml_select_next(1);
  }
  return (-1);
}

int axSelectPrev(){

  if(dataFormat==1){
    return xml_select_next(-1);
  }
  return (-1);
}

int axDeselect(){

  if(dataFormat==1){
    return xml_deselect();
  }
  return (-1);
}

int axRefine(char* expression){

  if(dataFormat==1){
    return xml_refine((xmlChar*)expression);
  }
  return (-1);
}

int axNamespace(char* namespacefURI){
  
  if ( axCheckState( AX_INIT, AX_OFF ) ) {
#ifdef DEBUG
    fprintf(stderr, "Error: parser in wrong state: axStateRec = %d\n", axStateRec);
#endif
    return -1;
  }
  if ( xmlStrlen( (const xmlChar*) namespacefURI ) > MAX_NS_LEN )
    {
#ifdef DEBUG
      fprintf(stderr, "Error: namespace URI greater than MAX_NS_LEN" );
#endif
      return -1;
    }
  else
    {
      strcpy ( (char*) namespaceURI, namespacefURI );
    }

  return(0);
}

int axLen(){

  if ( !gdata ) return(0);
  if ( gdata->charData ) return xmlStrlen( gdata->charData );

  return(0);
}

int axBuffer(int size){
  if ( size > 0 ) charbuffer = size;
  return(0);
}

char* axValue(){

  if ( !gdata ) return (char*)'\0';
  if ( !gdata->charData ) return (char*)'\0';

  return (char*)gdata->charData;
}

char* axBufferValue(){
  
  if ( !gdata ) return (char*)'\0';
  if ( !gdata->bufferPtr ) return (char*)'\0';
  
  if(((signed)xmlStrlen(gdata->bufferPtr) <= charbuffer) || (charbuffer == -1)){
    gdata->bufferData = gdata->bufferPtr;
    gdata->bufferPtr = NULL;
    return (char*)gdata->bufferData;
  }
  
  gdata->bufferData = (xmlChar*)xmlStrndup(gdata->bufferPtr, charbuffer);
  gdata->bufferPtr = (gdata->bufferPtr + charbuffer);
  
  return (char*)gdata->bufferData;
  
}

int axCache(int cachesize){
  maxEvalCacheSize = cachesize;
  return 0;
}

int control_load(int (*A)(), char *s){

  int retval = 0;
  char *p = NULL;

  do {
    p = s + xmlStrlen( (xmlChar*)s ) - 1;
    for ( ; p != s && ( *(p - 1) != ' ' || *(p - 1) != ',' ) ; p -- );
    if ( ((*A)( p )) ) {
      fprintf ( stderr, "Error loading %s\n", p );
      retval = -1;
    }
    for ( ; p != s && ( *(p - 1) != ' ' || *(p - 1) != ',' ) ; p -- );
    if ( p > s ) ( *p = '\0' );
  } while (p > s);

 return retval;
}

/* parse the AgentX control file (called iniFname) */

int axControl(char *iniFname, char *xmlFname )
{

    int retval = 0;
    dictionary * ini = NULL ;
    char *s = NULL;
    char *liniFname = NULL;

    /* if INI file not specified, use default */
    if ( ! iniFname || ! *iniFname ) {
        liniFname=(char *)malloc( MAX_CONFIG_FNAME_LENGTH ) ; 
        if( getINIFname( liniFname, MAX_CONFIG_FNAME_LENGTH ) < 0 ) ERROR
    }
    else
        liniFname = (char*) xmlStrdup ( (xmlChar*)iniFname );
  
#ifdef DEBUG
  printf("AgentX: INI file: %s\n", liniFname );
#endif

    /* load INI file */
    if( (ini = iniparser_load( liniFname )) == NULL ) {
        fprintf(stderr, "Error parsing %s\n", liniFname); ERROR
    }

    /* If parser not started, start the parser */
    if ( axCheckState( AX_INIT, AX_OFF ) ){
      if( axParserStart() < 0 ) {
	    fprintf(stderr, "Error starting parser\n"); ERROR
      }
    }

    /* set base URI, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:baseURI")) ) {

#ifdef DEBUG
  printf("AgentX: baseURI: %s\n", s );
#endif

        if ( axBaseUri(s) ) {
            fprintf(stderr, "Error setting base URI\n"); ERROR
        }
    }

    /* set format, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:format")) ) {
  
#ifdef DEBUG
  printf("AgentX: format: %s\n", s );
#endif

        if ( axFormat(atoi(s)) ) {
            fprintf(stderr, "Error setting format\n"); ERROR
        }
    }

    /* set cache size, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:cache")) ) {

#ifdef DEBUG
  printf("AgentX: cache: %s\n", s );
#endif

        if ( axCache(atoi(s)) ) {
	        fprintf(stderr, "Error setting cache size\n"); ERROR
        }
    }

    /* load ontology, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:ontology")) ) { 
  
#ifdef DEBUG
  printf("AgentX: ontology: %s\n", s );
#endif

        /* record if load successful, fail if not */
        if( (retval = control_load(axGetUri,s)) ) {
	        fprintf(stderr, "Error loading ontology document\n"); ERROR
        }
    }
  
    /* load mappings, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:map")) ) {

#ifdef DEBUG
  printf("AgentX: map: %s\n", s );
#endif

        /* record if load successful, fail if not */
        if( (retval = control_load(axGetUri,s)) ) {
	        fprintf(stderr, "Error loading mappings document\n"); ERROR
        }
    }
  
    /* load link, if specified */
    if( (s = iniparser_getstr(ini, "AgentX:link")) ) { 

#ifdef DEBUG
  printf("AgentX: link: %s\n", s );
#endif
        if( (retval = control_load(axGetUri,s)) ) {
	        fprintf(stderr, "Error loading link document\n"); ERROR
        }
    }  
  
    /* If no XML file specified, use default from INI file */
    if( xmlFname && *xmlFname ) 
        s = xmlFname;
    else 
        s = iniparser_getstr(ini, "AgentX:data");

    /* load data - data doc can also pull in ontology and mappings files... */
    if( s ) {

#ifdef DEBUG
  printf("AgentX: data: %s\n", s );
#endif
        if( (retval = control_load(axDataGetUri,s)) ) {
	        fprintf(stderr, "Error loading data document\n"); ERROR
        }
    }


 cleanup:
  iniparser_freedict(ini);
  if ( liniFname ) free ( liniFname );

  return retval ;
}


static int getINIFname( char *fname, size_t len)
{
    char *cp;

    if( (cp = getenv("HOME")) ) {
        /* $HOME/<CONFIG_DIRECTORY>/<CONFIG_FILENAME> */
        if( len < ( strlen(cp) + strlen(CONFIG_DIRECTORY) + strlen(CONFIG_FILENAME) + 3 ) ) {
            fprintf( stderr, "Home string is too long: %s\n", cp );
            return -1;
        }
        snprintf( fname, len, "%s/%s/%s", cp, CONFIG_DIRECTORY, CONFIG_FILENAME );
    }
    else {
        snprintf( fname, len, "agentx.ini" );
    }

    return 0;
}


int createWrap (struct xpathObjWrap **objWrap){

  int retval = 0;

  if (!objWrap) ERROR;

  if (!*objWrap) {
    *objWrap = (struct xpathObjWrap*)malloc(sizeof(struct xpathObjWrap));
    (*objWrap)->refs = 0;
    (*objWrap)->xpathObj = NULL;
    (*objWrap)->wrapTab = NULL;
    (*objWrap)->tabNr = 0;
  }

 cleanup:

  return (0);
}

int wrapObject (struct xpathObjWrap** objWrap, xmlXPathObjectPtr xpathObj){
 
  int retval = 0;
 
  if (!xpathObj) END;
  if (!objWrap) ERROR;

  if (!*objWrap) createWrap (objWrap);

  if ( (*objWrap)->xpathObj ) xmlXPathFreeObject( (*objWrap)->xpathObj );

  (*objWrap)->xpathObj = xpathObj;

 cleanup:
  
  return retval;
}

int incWrapRefs (struct xpathObjWrap *objWrap, int inc){

  int retval = 0, i = 0;

  if(!objWrap) ERROR;

  for (i = 0; i < objWrap->tabNr; i ++){
    if ( objWrap->wrapTab[i] )
      incWrapRefs( objWrap->wrapTab[i], inc );
  }

  objWrap->refs += inc;

 cleanup:

  return retval;
}

int linkWrap (struct xpathObjWrap **objWrap, struct xpathObjWrap *objWrapToAdd){

  int retval = 0;

  if ( !objWrapToAdd ) END;
  if (!objWrap) ERROR;

  if ( !*objWrap ) createWrap (objWrap);

  incWrapRefs( objWrapToAdd, 1 );
  (*objWrap)->tabNr ++;
  
  (*objWrap)->wrapTab = (struct xpathObjWrap**) realloc((*objWrap)->wrapTab,((*objWrap)->tabNr + 1) * sizeof(struct xpathObjWrap*));
  (*objWrap)->wrapTab[(*objWrap)->tabNr - 1] = objWrapToAdd;
  (*objWrap)->wrapTab[(*objWrap)->tabNr] = NULL;
  
 cleanup:
  
  return retval;
}

int axFreeWrap(struct xpathObjWrap** objWrap){

  int retval = 0, i = 0;

  if ( ! objWrap ) ERROR;
  if ( ! *objWrap ) END;

  for (i = 0; i < (*objWrap)->tabNr; i ++){
    axFreeWrap( &((*objWrap)->wrapTab[i]) );
  }

  if ( (*objWrap)->refs < 0 ){
    if ( (*objWrap)->xpathObj ) xmlXPathFreeObject( (*objWrap)->xpathObj );
    if ( (*objWrap)->wrapTab ) free((*objWrap)->wrapTab);
    free( *objWrap );
    *objWrap = NULL;
  }

 cleanup:
  
  return retval;
}

int freeWrap(struct xpathObjWrap** objWrap){

  int retval = 0;

  if ( ! objWrap ) ERROR;
  if ( ! *objWrap ) END;

  incWrapRefs( *objWrap, -1 );

  retval = axFreeWrap( objWrap );

 cleanup:

  return retval;

}

int axCountNodes (struct xpathObjWrap *objWrap, int *counter){

  int noLoc = 0, i = 0;

  if ( objWrap->xpathObj ){

    switch (objWrap->xpathObj->type){ 
      
    case(1):
      noLoc = objWrap->xpathObj->nodesetval->nodeNr;
      break;
      
    case(7):
      noLoc = ((xmlLocationSet*)objWrap->xpathObj->user)->locNr;
      break;

    default:
      noLoc = 1;
      break;
      
    }

    *counter += noLoc;
  }

  for ( i = 0; i < objWrap->tabNr; i++ ){
    axCountNodes (objWrap->wrapTab[i], counter);
  }

  return *counter;
}

int countNodes (struct xpathObjWrap *objWrap){
  
  int *counter = NULL, retval = 0;

  counter = (int*) malloc(sizeof(int));
  *counter = 0;
  
  retval = axCountNodes (objWrap, counter);

  free(counter);

  return retval;
}

int axLocateObject(struct xpathObjWrap *objWrap, int *no, xmlXPathObjectPtr *xpathObj){

  int retval = 0, i = 0;

  if(objWrap->xpathObj){

    if(*no == 1) {
      *xpathObj = objWrap->xpathObj;
    }
    
    *no --;
  }
  
  for (i = 0; i < objWrap->tabNr && *no > 0; i++){
    axLocateObject(objWrap->wrapTab[i], no, xpathObj);
  }
  
  return retval;
}

int locateObject(struct xpathObjWrap *objWrap, int no, xmlXPathObjectPtr *xpathObj){

  int retval = 0;
  int *startNo = NULL;

  if (no < 1) ERROR;

  startNo = (int*) malloc (sizeof(int));
  *startNo = no;
  *xpathObj = NULL;

  retval = axLocateObject(objWrap, startNo, xpathObj);

  free(startNo);

 cleanup:

  return retval;
}

int axLocateNode(struct xpathObjWrap *objWrap, int *no, xmlNodePtr *node, int *start, int *end){

  int noLoc = 0, retval = 0, i = 0;
  xmlLocationSet* xmlLS = NULL;

  if (!objWrap) ERROR;
  if (!node) ERROR;
  if (*no < 1) ERROR;
  
  if(objWrap->xpathObj){

    switch (objWrap->xpathObj->type){ 
      
    case(1):
      noLoc = objWrap->xpathObj->nodesetval->nodeNr;
      if( *no <= noLoc){
	*node = objWrap->xpathObj->nodesetval->nodeTab[*no - 1];
	if (start) *start = -1;
	if (end) *end = -1;
      }
      break;
      
    case(7):
      noLoc = ((xmlLocationSet*)objWrap->xpathObj->user)->locNr;
      if ( *no <= noLoc ){
	xmlLS = (xmlLocationSet*)objWrap->xpathObj->user;
	*node = xmlLS->locTab[ *no - 1]->user;
	if (start) *start = xmlLS->locTab[ *no - 1]->index;
	if (end) *end = xmlLS->locTab[ *no - 1]->index2;
      }
      break;

    default:
      noLoc = 1;
      if( *no <= noLoc){
	*node = objWrap->xpathObj->user;
	if (start) *start = objWrap->xpathObj->index;
	if (end) *end = objWrap->xpathObj->index2;
      }
      break;

    }

    (*no) -= noLoc;

  }

  for  ( i = 0; i < objWrap->tabNr && *no > 0; i++ ){
    if ( axLocateNode(objWrap->wrapTab[i], no, node, start, end) ) ERROR;
  }
  
 cleanup:
  
  return retval;
}

int locateNode(struct xpathObjWrap *objWrap, int no, xmlNodePtr *node, int *start, int *end){

  int retval = 0;
  int *startNo = NULL;

  startNo = (int*) malloc (sizeof(int));
  *startNo = no;
  *node = NULL;

  if (start) *start = -1;
  if (end) *end = -1;

  retval = axLocateNode(objWrap, startNo, node, start, end);

  free(startNo);

  return retval;
}
  
int axCurrent(){

  struct xml_data* xml_data_local=NULL;
  struct ctxStr* ctxStrLocal=NULL;
  int total=0;
  int current=0;
  int retval=0;

  if(!ctxStrPtr)return(-1);
  if(!gdata)return(-1);

  xml_data_local=ctxStrPtr->xml_data_origin;
  ctxStrLocal=ctxStrOrigin;

  while(xml_data_local){
    if(xml_data_local->prev){
      printf("%s %d %d\n",(char*)xml_data_local->class,xml_data_local->data_sets_position,xml_data_local->data_sets_count);
      if(xml_data_local->set==&xml_data_local->item_set)printf("%s %d %d\n",(char*)xml_data_local->property,xml_data_local->item_pos,xml_data_local->item_count);
    }
    else{
      while(ctxStrLocal){
	if(ctxStrLocal->xml_data_origin){
	  if(ldata==gdata){
	    current=total+ctxStrLocal->xml_data_origin->data_sets_position;
	  }
	  total+=ctxStrLocal->xml_data_origin->data_sets_count;
	}
	ctxStrLocal=ctxStrLocal->next;
      }
      printf("%s %d %d\n",(char*)ctxStrPtr->xml_data_origin->class,current,total);
      if(xml_data_local->set==&xml_data_local->item_set)printf("%s %d %d\n",(char*)xml_data_local->property,xml_data_local->item_pos,xml_data_local->item_count);
    }
    xml_data_local=xml_data_local->next;
  }
  
  return(retval);
}

int registerStandardDelimiters(xmlChar*** delimiters){

  int noStdDels = 1;
  int retval = 0, size = 0;
  
  if ( !delimiters ) ERROR;

  if (*delimiters) for (size = 0; (*delimiters)[size]; size++);

  *delimiters = (xmlChar**) realloc( *delimiters, (size + noStdDels + 1) * sizeof(xmlChar*) );
  (*delimiters)[size + noStdDels] = NULL;

  (*delimiters)[size + noStdDels - 1] = xmlStrdup ( (xmlChar*)"\n" );
  
 cleanup:
  
  return retval;
}

int resolveExpr(xmlChar* expr,xmlChar** resExpr){

  xmlChar* pos=NULL;
  int i=0,j=0,k=0;
  int size=0;
  int retval=0;
  int expressionSize=0;
  xmlChar* expressionStart=NULL;
  xmlChar* expression=NULL;
  float evaluation=0;
  float feval=0;

  /* replace ||string|| with resolveExpression(string)
   */

  if(!expr){return(-1);}
  if(!resExpr){return(-1);}

  if(*resExpr)free(*resExpr);
  *resExpr=(xmlChar*)xmlStrdup(expr);

  if(!xmlStrstr(*resExpr,(const xmlChar*)"||")){return(0);}

  /* read through expr until null character is reached */

  for(i=0;*(expr+i)!='\0';i++){

    if((*(expr+i)=='|')&&(*(expr+i+1)=='|')){
      
      expressionStart=expr+i+2;
      expressionSize=(int)(xmlStrstr((const xmlChar*)(expr+i+2),(const xmlChar*)"||")-(expr+i+2));
      if(expression)free(expression);
      expression=(xmlChar*)malloc((expressionSize+2)*sizeof(xmlChar));

      /* expressionSize must be less than 50 */

      if(expressionSize>=50) ERROR;

      xmlStrPrintf((xmlChar*)expression,(expressionSize+1),(const xmlChar*)"%s",(char*)expressionStart);

      evaluation=evaluateExpression((char*)expression);

      feval = (float)((int) evaluation);
      if (feval != evaluation) {
	if (evaluation > 0)
	  evaluation = feval + 1;
	else {
	  if (evaluation < 0 && feval == 0)
	    evaluation = 0;
	  else
	    evaluation = feval;
	}
      }

      if(pos)free(pos);
      pos=(xmlChar*)malloc(100*sizeof(xmlChar));

      sprintf((char*)pos,"%d",(int)evaluation);

      if(pos){
	for(j=0;*(pos+j)!='\0';j++){
	  if(k>=size-2){
	    *resExpr=(xmlChar*)realloc(*resExpr,(size+100)*sizeof(xmlChar));
	    size+=100;
	  }
	  *(*resExpr+k)=*(pos+j);
	  k++;
	}
      }
      
      i+=expressionSize+3;
    }
    
    else{
      *(*resExpr+k)=*(expr+i);
      k++;
    }
  }
  
  *(*resExpr+k)='\0';

 cleanup:

  if(pos)free(pos);
  if(expression)free(expression);
  if(retval==-1){
    if(*resExpr)free(*resExpr);
  }

  return retval;
}

int axSelectNo(int pos){

  int no = 0, retval = 0;

  if ( !gdata ) ERROR;
  if ( pos < 1 ) ERROR;

  if (gdata->set == &gdata->data_sets){
    if (pos > gdata->data_sets_count) ERROR;
    gdata->data_sets_position = pos;
    locateNode(gdata->data_sets, gdata->data_sets_position, &gdata->sel, &gdata->sel_index, &gdata->sel_index2);
    return (gdata->data_sets_count - pos + 1);
  }

  else{
    if (pos > gdata->item_count) ERROR;
    no = pos - gdata->item_pos;
    if (no < 0) {
      while (no){
	axSelectPrev();
	no ++;
      }
    }
    else{
      while (no){
	axSelectNext();
	no --;
      }
    }

    return (gdata->item_count - pos + 1);
  }
  
 cleanup:
  
  return retval;
}

char* axAbout(char* name){

  int ssize=0,nods=0,noProp=0;
  char *retval='\0',*rstr=NULL;

  if(!gdata)return retval;

  if((nods=axSelectType("Metadata",1))<1){
    return retval;
  }

  ssize=xmlStrlen((const xmlChar*)name)+6;
  rstr=(char*)malloc(ssize*sizeof(char));
  xmlStrPrintf((xmlChar*)rstr,ssize,(const xmlChar*)"name=%s",name);
  
  if(axRefine(rstr)<1){
    axDeselect();
    goto cleanup;
  }

  if((noProp=xml_get((xmlChar*)"content"))<1)goto cleanup;

  retval=axValue();
  axDeselect();
  axDeselect();

 cleanup:
  
  if(rstr)free(rstr);

  return retval;
}

char* axSelectValue(char* property){

  int noProp = 0;
  char *value = NULL;

  if( ( noProp = xml_get( (xmlChar*)property ) ) > 0 ){
    value = axValue();
    axDeselect();
    return value;
  }
  
  return (char*)'\0';
}

/* ax XPath functions */


void xmlBaseFunction(xmlXPathParserContextPtr ctxt, int nargs){

  xmlNodePtr node = NULL;
  xmlXPathObjectPtr strobj = NULL;

  node=ctxt->context->node;
  strobj = xmlXPathWrapString(xmlNodeGetBase(node->doc,node));
  valuePush(ctxt,strobj);

}


int registerInternalFunctions(xmlXPathContextPtr ctxt){

  xmlXPathRegisterFunc(ctxt, (const xmlChar *)"xmlbase",
		       xmlBaseFunction);
  return(0);
}


int resolveVars(struct ctxStr* ctxStrLocal, xmlChar* expr, xmlChar** resExpr){

  /* replace variables in expr (indicated by a $) to produce a valid XPointer expression */

  int i = 0, j = 0, k = 0;
  xmlChar* pos = NULL;
  int size = 0;
  int retval = 0;

  if (!expr) return(-1);
  if (!ctxStrLocal) return(-1);
  if (!resExpr) return(-1);

  if (*resExpr) free(*resExpr);
  *resExpr = (xmlChar*)xmlStrdup(expr);

  if(!xmlStrchr(expr, '$')){return(0);}

  /* read through expr until null character is reached */

  for(i = 0; *(expr + i) != '\0'; i++){

    /*
       replace $1 with ctxStrLocal->var
       replace $2 with ctxStrLocal->var2
       replace $3 with ldata->data_sets_position
    */
    
    if(*(expr + i) == '$'){
      
      if (pos) free(pos);
      pos = NULL;
      
      switch(*(expr + i + 1)){
	
      case('1'):
	if(ctxStrLocal->var)pos=xmlStrdup(ctxStrLocal->var);
	break;
	
      case('2'):
	if(ctxStrLocal->var2)pos=xmlStrdup(ctxStrLocal->var2);
	break;
	
      case('3'):
	if(ldata){
	  pos=(xmlChar*)malloc(100*sizeof(xmlChar));
	  sprintf((char*)pos,"%d",*(ldata->set_pos));
	}
	break;

      default:
	pos = (xmlChar*) xmlStrndup(expr+i,2);

      }

      for(j = 0; *(pos + j) != '\0'; j++){
	if(k >= size - 2){
	  *resExpr = (xmlChar*)realloc(*resExpr, (size+100) * sizeof(xmlChar));
	  size += 100;
	}
	*(*resExpr + k) = *(pos + j);
 	k++;
      }
          
      i += 1;
    }
    
    else{

      if(k >= size - 2){
	*resExpr = (xmlChar*)realloc(*resExpr, (size + 100) * sizeof(xmlChar));
	size += 100;
      }
      
      *(*resExpr + k) = *(expr + i);
      k++;
    }

  }

  *(*resExpr + k) = '\0';
  if (pos) free(pos);

  return (retval);
}


int toHash(xmlChar* str){

  int hash = 0;

  if (!str) return (0);

  while (*str) {
    hash = hash + (2 * (*str));
    str ++;
  }

  return hash;
}
    

int xptr_query(struct ctxStr* ctxStrLocal, xmlNodePtr context, xmlChar* xptrExpr){   /* evaluate an XPointer query */

  void* load = NULL;
  xmlChar* resExpr = NULL;
  int retval = 0, expHash = 0, hash = 0;
  struct evalCache *evalLocal = NULL;
  struct cacheBtree *treeLocal = NULL, *treeLocalPrev = NULL;
  xmlXPathObjectPtr object = NULL;
  xmlNodePtr root_prev = NULL, root_next = NULL, root_parent = NULL;
  xmlChar* evalContext = NULL;
  struct xpathObjWrap *objWrap = NULL;
  int found = 0;

  if(!ctxStrLocal){return(-1);}
  if(!ctxStrLocal->xpathCtx){return(-1);}
  if(!context){return(-1);}
  if(context->type!=XML_ELEMENT_NODE){return(-1);}
  if(!xptrExpr){return(-1);}

  if (resolveVars (ctxStrLocal,xptrExpr,&resExpr)) {return(-1);}

  if(linksResolved){
    linkup(context);
    getevalcontext(context,&evalContext);
  }

  else{
    evalContext=(xmlChar*)malloc(20*sizeof(xmlChar));
    sprintf((char*)evalContext,"%d",(int)context);
  }

  /* check cache for evaluation */
  
  expHash = toHash(resExpr);

  if (useCache == 1){
    
    searchCache(xptrCacheOrigin, &treeLocalPrev, resExpr, evalContext, &load, &found);

      if ( found ) {
	freeWrap( &(ctxStrLocal->xpathObj) );
	ctxStrLocal->xpathObj = (struct xpathObjWrap*) load;
	incWrapRefs(ctxStrLocal->xpathObj, 1);
	if (linksResolved) axunlink(context);
	free (evalContext);
	free (resExpr);
	return (0);
      }
  }

  /* set document root node */
  
  ctxStrLocal->xpathCtx->doc->children=context;
  root_prev=context->prev;
  context->prev=NULL;
  root_next=context->next;
  context->next=NULL;
  root_parent=context->parent;
  context->parent=NULL;

  object = xmlXPtrEval(resExpr,ctxStrLocal->xpathCtx);

  /* make sure there is an xmlXPathObject */

  if (!object) END

  /* check the returned object to make sure it contains some useful data */

  else if((object->type==XPATH_NODESET)&&(!object->nodesetval)){
    xmlXPathFreeObject(object);
    object=NULL;
  }

  else if((object->type==XPATH_NODESET)&&(object->nodesetval->nodeNr==0)){
    xmlXPathFreeObject(object);
    object=NULL;
  }

  else if((object->type==XPATH_NODESET)&&(!object->nodesetval->nodeTab)){
    xmlXPathFreeObject(object);
    object=NULL;
  }

  else if((object->type==XPATH_LOCATIONSET)&&(!object->user)){
    xmlXPathFreeObject(object);
    object=NULL;
  }

 cleanup:

  /* restore document root node */
  
  ctxStrLocal->xpathCtx->doc->children->parent=root_parent;
  ctxStrLocal->xpathCtx->doc->children->prev=root_prev;
  ctxStrLocal->xpathCtx->doc->children->next=root_next;
  ctxStrLocal->xpathCtx->doc->children=ctxStrLocal->original_doc_root;

  wrapObject(&objWrap, object);

  /* add evaluation to cache */

  if (useCache == 1){

    addToCache (&treeLocalPrev, resExpr, evalContext, objWrap);
    if (!xptrCacheOrigin) xptrCacheOrigin = treeLocalPrev;
    if (objWrap) objWrap->refs++;
    free(resExpr);
    free(evalContext);

  }

  freeWrap(&(ctxStrLocal->xpathObj));
  ctxStrLocal->xpathObj = objWrap;
  incWrapRefs(objWrap, 1);

  if(linksResolved){
    axunlink(context);
  }

  return (retval);

}

/* select next (dir > 0), select prev (dir < 0)
 */

int xml_select_next(int dir){

  int retval = 0, pos = 0;
  struct ctxStr* ctxStrLocal = NULL;

  if ( dir < 0 ) dir = -1;
  else dir = 1;
  
  ctxStrLocal = ctxStrPtr;

  while(ctxStrLocal){

    if(ldata){
      
      if(ldata->set == &(ldata->item_set)){
	
	retval = xml_get_next(dir);
	
	return retval;
      }
      
      else{

	pos = ldata->data_sets_position + dir;
	if((pos <= ldata->data_sets_count) && (pos > 0)){

	  ldata->data_sets_position = pos;

	  if(!ldata->data_sets_type){
	    locateNode(ldata->data_sets, ldata->data_sets_position, &ldata->sel, &ldata->sel_index, &ldata->sel_index2);
	  }

	  ctxStrPtr = ctxStrLocal;
	  return(gdata->data_sets_count - gdata->data_sets_position + 1);
	}
	else{

	  if(ldata->prev){
	    return(0);
	  }

	}
      }
    }

    if (dir < 0) ctxStrLocal = ctxStrLocal->prev;
    else ctxStrLocal = ctxStrLocal->next;

  }

  return(0);
}


int xml_refine(const xmlChar* expression){

  int retval = 0;
  xmlChar* operand1 = NULL;
  xmlChar* operand2 = NULL;
  xmlChar* orig_operator = NULL;
  xmlChar* operator = NULL;
  int lenop = 1, i = 0, j = 0, k = 0;
  xmlXPathObjectPtr newObject = NULL;
  struct ctxStr* ctxStrLocal = NULL;
  xmlNodeSetPtr nodeSet = NULL;

  
  if (!expression) return(-1);
  if (!ctxStrPtr) return(-1);
  if (!gdata) return(-1);
  if (!*(gdata->set)) return(-1);

  /* look for operator in expression */

  if (!orig_operator) orig_operator 
			= (xmlChar*) xmlStrstr(expression, (const xmlChar*) "<" );
  if (!orig_operator) orig_operator 
			= (xmlChar*) xmlStrstr(expression, (const xmlChar*) "=" );
  if (!orig_operator) orig_operator 
			= (xmlChar*) xmlStrstr(expression, (const xmlChar*) ">" );

  if (!orig_operator) {
    orig_operator = (xmlChar*) xmlStrstr(expression, (const xmlChar*) "<=" );
    lenop = 2;
  }

  if (!orig_operator) {
    orig_operator = (xmlChar*) xmlStrstr(expression, (const xmlChar*) ">=" );
    lenop = 2;
  }
  
  if (!orig_operator) return(-1);

  /* extract operator */

  operator = (xmlChar*) xmlStrndup(orig_operator, lenop);

  /* extract operands */

  operand1 = (xmlChar*) xmlStrndup(expression, (int)(orig_operator - expression));

  operand2 = (xmlChar*)(orig_operator + lenop);

  if ((!operand1) || (!operand2)) ERROR;

  ctxStrLocal = ctxStrPtr;

  nodeSet = (xmlNodeSetPtr) malloc(sizeof(xmlNodeSet));
  nodeSet->nodeTab = NULL;

  /* the nodes of interest */

  i = 1;
  locateNode(*(ldata->set), i, &ldata->sel, &ldata->sel_index, &ldata->sel_index2);
  k = *(ldata->set_pos);
  *(ldata->set_pos) = i;

  for(i = 2; ldata->sel; i++){

    /* deal with the = operator */

    if ( ! xmlStrcmp(operator,(const xmlChar*)"=") ){
      if ( xml_get(operand1) > 0 ){
	if( !xmlStrcmp((xmlChar*)axValue(), operand2) ){
	  xml_deselect();
	  j++;
	  nodeSet->nodeTab = realloc(nodeSet->nodeTab, (j+1) * sizeof(xmlNodePtr));
	  *(nodeSet->nodeTab + j - 1) = ldata->sel;
	  *(nodeSet->nodeTab + j) = NULL;
	  nodeSet->nodeNr = j;
	}
	else{
	  xml_deselect();
	}
      }
    }

    /* deal with other operators */
    /* need to complete this */

    locateNode(*(ldata->set), i, &ldata->sel, &ldata->sel_index, &ldata->sel_index2);
    *(ldata->set_pos) = i;

  }

  locateNode(*(ldata->set), k, &ldata->sel, &ldata->sel_index, &ldata->sel_index2);
  *(ldata->set_pos) = k;

  /* create a new object and replace ldata->data_sets with that object */

  if(nodeSet->nodeTab){

    create_xml_data(ctxStrLocal);

    newObject = (xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));
    newObject->nodesetval = nodeSet;
    newObject->type = 1;
    newObject->boolval = 0;

    freeWrap( &ldata->data_sets );
    ldata->data_sets = NULL;
    wrapObject( &ldata->data_sets, newObject );

    ldata->class = xmlStrdup( ldata->prev->class );
    ldata->data_sets_position = 1;
    ldata->data_sets_count = j;
    ldata->set = &ldata->data_sets;
    ldata->set_pos = &ldata->data_sets_position;
    locateNode(ldata->data_sets, ldata->data_sets_position, &ldata->sel, &ldata->sel_index, &ldata->sel_index2);

  }

  else{
    if (nodeSet) xmlXPathFreeNodeSet(nodeSet);
  }

 cleanup:

  if (operand1) free(operand1);
  if (operator) free(operator);
  if (retval == -1) j = retval;

  return j;

}

int xml_deselect(){
 
  struct xml_data* xml_data_local=NULL;
  int retval=0;

  if(!ctxStrPtr){return(-1);}
  if(!gdata){return(-1);}

  if(gdata->set==&gdata->item_set){
    gdata->set=&gdata->data_sets;
    gdata->set_pos=&gdata->data_sets_position;
    locateNode(gdata->data_sets,gdata->data_sets_position,&gdata->sel,&gdata->sel_index,&gdata->sel_index2);
    return (0);
  }

  if(!gdata->prev){

    ctxStrPtr=ctxStrOrigin;

    while(ctxStrPtr){

      if(gdata){

	if(destroy_xml_data(ctxStrPtr)){return(-1);}

	free(gdata);
	gdata=NULL;
	ctxStrPtr->xml_data_origin=NULL;
      }
      ctxStrPtr=ctxStrPtr->next;
    }
    ctxStrPtr=ctxStrOrigin;
  }
  
  else{

    if(destroy_xml_data(ctxStrPtr)){return(-1);}
    
    xml_data_local=gdata->prev;
    free(gdata);
    gdata=xml_data_local;
    gdata->next=NULL;
  }

  return(retval);
}

int create_xml_data(struct ctxStr* ctxStrLocal){

  int retval=0;

  if (!ctxStrLocal) return(-1);

  if (!ldata){
    ldata = ctxStrLocal->xml_data_origin = (struct xml_data*) malloc(sizeof(struct xml_data));
    ldata->prev = NULL;
  }
  
  else{
    ldata->next = (struct xml_data*) malloc(sizeof(struct xml_data));
    ldata->next->prev = ldata;
    ldata = ldata->next;
  }
  
  ldata->next = NULL;
  ldata->property = NULL;
  ldata->xlocator = NULL;
  
  ldata->cdata_list = NULL;
  ldata->data_sets = NULL;
  ldata->data_sets_position = 1;
  ldata->data_sets_count = 0;
  ldata->data_sets_type = 0;
  ldata->item_set = NULL;
  ldata->item_set_position = 1;
  ldata->item_pos = 1;
  ldata->item_count = 0;
  ldata->startIndex = NULL;
  ldata->endIndex = NULL;
  ldata->skip_delimiter = 1;
  ldata->index = 1;
  ldata->set = &ldata->data_sets;
  ldata->set_pos = &ldata->data_sets_position;

  ldata->delimiter = NULL;
  ldata->charPtr = NULL;
  ldata->charData = NULL;
  ldata->freeable = NULL;
  ldata->class = NULL;
  ldata->next_resource = NULL;
  ldata->bufferPtr = NULL;
  ldata->bufferData = NULL;
  ldata->sel = NULL;
  ldata->sel_index = -1;
  ldata->sel_index2 = -1;
  ldata->card = -1;
  ldata->data_sets_type = 0;
  ldata->startChar = 1;
  ldata->endChar = -1;

  return(retval);
  
}

int destroy_xml_data(struct ctxStr* ctxStrLocal){

  int i=0,retval=0;

  if(!ctxStrLocal)return(-1);
  if(!ldata)return(-1);
  if(ldata->class)free(ldata->class);
  if(ldata->next_resource)free(ldata->next_resource);
  if(ldata->data_sets)freeWrap(&ldata->data_sets);
  if(ldata->item_set)freeWrap(&ldata->item_set);
  if(ldata->property)free(ldata->property);
  
  if(ldata->cdata_list){
    if(ldata->freeable){
      for(i=0;*(ldata->freeable+i);i++){
	free(*(ldata->freeable+i));
      }
    }
    free(ldata->cdata_list);
    ldata->cdata_list=NULL;
    ldata->freeable=NULL;
  }

  if ( charData == &(ldata->charData) ) { free (*charData); charData = NULL; }

  return (retval);
 
}


int axParserFinish(){

  struct ctxStr* ctxStrLocal=NULL;
  struct xml_data* xml_data_local=NULL;

  ctxStrPtr=ctxStrOrigin;

  /* free ctxStr structs */

  while(ctxStrPtr){

    if (ctxStrPtr->xpathObj) freeWrap(&ctxStrPtr->xpathObj);

    /* free docs and contexts */
    
    if (ctxStrPtr->next_resource) free(ctxStrPtr->next_resource);
    if (ctxStrPtr->var) free(ctxStrPtr->var);
    if (ctxStrPtr->var2) free(ctxStrPtr->var2);
    if (ctxStrPtr->copyNodes) freeWrap (&ctxStrPtr->copyNodes);

    gdata = ctxStrPtr->xml_data_origin;

    /* free xml_data_structures */

    while(gdata){

      destroy_xml_data(ctxStrPtr);
      xml_data_local = gdata->next;
      free(gdata);
      gdata = xml_data_local;

    }

    if (ctxStrPtr->xpathCtx){

      if (ctxStrPtr->xpathCtx->doc) xmlFreeDoc(ctxStrPtr->xpathCtx->doc);
      xmlXPathFreeContext (ctxStrPtr->xpathCtx);

    }

    ctxStrLocal = ctxStrPtr->next;
    free(ctxStrPtr);
    ctxStrPtr = ctxStrLocal;
  }

  ctxStrOrigin = NULL;

  /* free static chars */

  if (base_uri_string){
    free(base_uri_string);
    base_uri_string = NULL;
  }
  if (namespaceURI) {
    free(namespaceURI);
    namespaceURI = NULL;
  }

  /* need to clear cache */


  /* reset parser internal state */

  axStateRec = 0;
  xmlCleanupParser();
  
  return(0);

}

xmlChar* locationCastToString(xmlXPathObjectPtr location){

  /* used for ranges and points */

  xmlChar* buffer=NULL;
  xmlNodePtr node=NULL,node2=NULL;
  int index=0,index2=0,i=0;

  if(!location)return (xmlChar*)'\0';

  node=(xmlNodePtr)location->user;
  node2=(xmlNodePtr)location->user2;
  index=location->index;
  index2=location->index2;

  if(location->type==XPATH_POINT){
    return (xmlChar*)'\0';
  }
  
  if(location->type==XPATH_RANGE){

    if(index<0)index=0;
    if(index2<0)index2=0;
    if((!node)||(!node2))return (xmlChar*)'\0';
    buffer=(xmlChar*)malloc((1+index2-index)*sizeof(char));

    /* node and node2 must be the same node */

    for(i=index;i<index2;i++){
      *(buffer+i-index)=*(node->children->content+i);
    }

    *(buffer+i-index)='\0';
    return buffer;
  }
  
  return (xmlChar*)'\0';
  
}

int axBaseUri(char* local_base_uri_string){
  
  if ( axCheckState( AX_INIT, AX_OFF ) ) {
#ifdef DEBUG
    fprintf(stderr, "Error: parser in wrong state: axStateRec = %d\n", axStateRec);
#endif
    return -1;
  }
  if ( xmlStrlen( (const xmlChar*) local_base_uri_string ) > MAX_NS_LEN ) 
    {
#ifdef DEBUG
      fprintf(stderr, "Error: base URI greater than MAX_NS_LEN" );
#endif
      return -1;
    }
  else
    {
      strcpy ( (char*) base_uri_string, local_base_uri_string );
    }

  return(0);
}

int type_query(int type,xmlNodePtr context,xmlChar* xptrExpr){

  /*
    behaviour type codes:
    0 = query any docs
    1 = query data docs
    2 = query non-data docs loaded from data docs
    3 = query non-data docs loaded from API
    4 = query last document
    5 = query current document
    6 = query non-data docs loaded from current doc
  */
  
  /* non-data docs are mapping, ontology or rules documents */

  int retval=0;
  xmlNodePtr evalctxt=NULL;
  struct ctxStr* ctxStrLocal2=NULL;

  evalctxt=context;
  ctxStrLocal2=ctxStrOrigin;

  while(ctxStrLocal2){

    if((ctxStrLocal2->type==type)||(type==0)||((type==5)&&(ctxStrLocal2==ctxStrPtr))
       ||((type==4)&&(!ctxStrLocal2->next))||((type==6)&&(ctxStrLocal2->parent==ctxStrPtr))){

      if(!context)evalctxt=ctxStrLocal2->original_doc_root;
      if(xptr_query(ctxStrLocal2,evalctxt,xptrExpr))return(-1);

    }

    else{
      if(ctxStrLocal2->xpathObj){
	freeWrap(&(ctxStrLocal2->xpathObj));
      }
      ctxStrLocal2->xpathObj=NULL;
    }

    ctxStrLocal2=ctxStrLocal2->next;
  }

  return retval;
}

int getDataElements(struct xpathObjWrap *objWrap, xmlChar*** elementList, int* elementCount, xmlChar*** freeable){
 
  int counter = 0;
  int i = 0, j = 0;
  xmlNodePtr node = NULL;
  xmlLocationSetPtr location_set = NULL;
  xmlXPathObjectPtr location = NULL;
  xmlXPathObjectPtr xpathObj = NULL;
  
  if ( !objWrap ) return(-1);
  if ( !elementList ) return(-1);

  locateObject(objWrap, 1, &xpathObj);

  for( i = 2; xpathObj; i++ ){
    
    /* deal with node sets */
    
    if( xpathObj->type == XPATH_NODESET ){
      
      for( j = 0; j < xpathObj->nodesetval->nodeNr; j++){
	
	node = xpathObj->nodesetval->nodeTab[j];
	
	counter = counter + 1;
	
	if (elementList){
	  
	  *elementList = realloc(*elementList, ((counter) + 1) * sizeof(xmlChar*));
	  
	  if ( node->children ) {
	    *((*elementList) + (counter) - 1) = node->children->content;
	  } 
	  else {
	    *((*elementList) + (counter) - 1) = NULL;
	  }
	  
	  *((*elementList) + counter) = NULL;
	}
      }
    }
    
    /* deal with location sets */

    if ( xpathObj->type == XPATH_LOCATIONSET){
      
      location_set = (xmlLocationSetPtr)xpathObj->user;
      
      for ( j = 0; j < location_set->locNr; j++){
	
	location = location_set->locTab[j];
	if(((xmlNodePtr)location->user)->children){
	  counter=counter+1;
	  
	  if(elementList){
	    *elementList=realloc(*elementList,((counter)+1)*sizeof(xmlChar*));
	    *((*elementList)+(counter)-1)=locationCastToString(location);
	    *((*elementList)+counter)=NULL;
	    
	  }
	}
      }
      
      if(freeable){
	*(freeable)=*elementList;
      }
      
    }

    locateObject(objWrap, i, &xpathObj);
  }

  if (elementCount) *elementCount = counter;

  return(0);
}


int collectResults (struct xpathObjWrap **results){

  struct ctxStr* ctxStrLocal=NULL;

  if (!results) return(-1);

  ctxStrLocal = ctxStrOrigin;

  while (ctxStrLocal){
    if (ctxStrLocal->xpathObj) linkWrap(results, ctxStrLocal->xpathObj);
    freeWrap(&ctxStrLocal->xpathObj);
    ctxStrLocal->xpathObj = NULL;
    ctxStrLocal = ctxStrLocal->next;
  }

  return 0;
}


int xml_get_next(int dir){

  int i = 0, retval = 0, counter = 0, imax = 0, pos = 0;
  xmlChar *begin = NULL, *end = NULL, **istr = NULL;
  xmlChar* cptr = NULL;

  if (!ctxStrPtr) return (0);
  if (!gdata) return (0);
  if (!gdata->cdata_list) return (0);

  if ( dir < 0 ) {
    dir = -1;
  }
  else{
    dir = 1;
  }
  counter = 2 * dir;
  
  pos = gdata->item_pos + dir;
  if ((pos <= 0 ) || (pos > gdata->item_count)) return (0);
  
  if (! gdata->delimiter){
    
    gdata->item_pos += dir;
    gdata->index += dir;
    gdata->charPtr = *( gdata->cdata_list + gdata->item_pos - 1);
    if ( charData ) { free (*charData); *charData = NULL; charData = NULL;}
    gdata->charData = gdata->charPtr;
    gdata->bufferPtr = gdata->charData;
    retval = gdata->item_count - gdata->item_pos + 1;
    gdata->item_set_position = gdata->item_pos;
    gdata->item_count = retval;
    
  }
  
  else{
    
    istr = gdata->cdata_list;
    imax = countNodes(gdata->item_set);
    cptr = gdata->charPtr;

    for(i = gdata->item_set_position - 1; i < imax && i >= 0 && counter; i += dir){
      
      ax_strstr(istr[i],
		cptr,
		gdata->delimiter,
		&counter,
		&begin,
		&end,
		gdata->startChar,
		gdata->endChar);

      if (i < imax && i >= 0) cptr = istr[i+dir];

    }

    if (counter) return(0);
    
    if (*(end + 1)){
      
      if ( charData ){ free (*charData); *charData = NULL;}
      gdata->charData = xmlStrndup(begin, ((end - begin) + 1));
      charData = &(gdata->charData);
      
      gdata->charPtr = begin;
      gdata->bufferPtr = gdata->charData;
      
    }
    
    else{
      
      gdata->charPtr = begin;
      if (charData) { free (*charData); *charData = NULL; charData = NULL; }
      gdata->charData = begin;
      gdata->bufferPtr = begin;
      
    }
    
    gdata->item_set_position = i - dir + 1;
    gdata->item_pos += dir;
    gdata->index += dir;
  }
  
  return (gdata->item_count - gdata->item_pos + 1);
  
}


int charMalloc(void **memptr,int allocSize,int reqSize){
  
  if(!memptr)return(-1);
  
  if(*memptr){
    if(allocSize<reqSize){
      free(*memptr);
      *memptr=malloc(reqSize*sizeof(xmlChar));
    }
  }
  else *memptr=malloc(reqSize*sizeof(xmlChar));

  return (0);
}


int ax_sprintf(xmlChar **memptr,int allocSize,int reqSize,char *format, ...){

  va_list args;

  charMalloc((void**)(memptr),allocSize,reqSize);
  va_start(args,format);
  vsprintf((char*)(*memptr),format,args);

  return(0);
}


int nextSiblingOfType(xmlNodePtr* node,unsigned int type){

  /* if **node is of the correct type
     keep, otherwise move to the next 
     sibling node of the correct type
  */

  if(!node)return(-1);

  if(*node){
    while((*node)->type!=type){
      *node=(*node)->next;
      if(!(*node))break;
    }
  }
  return(0);
}


int startResource(struct xpathObjWrap **objWrap, xmlChar *resource){

  int slen = 0, retval = 0;
  xmlChar* xptrExpr = NULL;

  if (!resource) return(-1);
  if (!objWrap) return(-1);

  xptrExpr = (xmlChar*)malloc(1000 * sizeof(xmlChar));

  slen = 128 + (3 * xmlStrlen(resource));
  ax_sprintf(&xptrExpr, 1000, slen, "xpointer(//*[@rdf:ID=substring-after('%s',concat(xmlbase(),'#')) or @rdf:about='%s' or @rdf:about=substring-after('%s',xmlbase())])", resource, resource, resource);

  if (type_query(5, NULL, xptrExpr)) ERROR;
  collectResults(objWrap);
  if (type_query(6, NULL, xptrExpr)) ERROR;
  collectResults(objWrap);
  if (type_query(3, NULL ,xptrExpr)) ERROR;
  collectResults(objWrap);

 cleanup:
  if (xptrExpr) free(xptrExpr);

  return(0);

}

int resourceURI(xmlChar** resource,xmlNodePtr node){

  xmlChar *nodeBase=NULL,*frag=NULL;
  int slen=0,retval=0;
  
  nodeBase=xmlNodeGetBase(node->doc,node);
  frag=getProperty(node,(xmlChar*)"about",(xmlChar*)"http://www.w3.org/1999/02/22-rdf-syntax-ns#");
  
  if(frag){
    if(*frag=='#'){
      slen=xmlStrlen(frag)+xmlStrlen(nodeBase)+1;
      ax_sprintf(resource,0,slen,"%s%s",nodeBase,frag);
    }
    else{
      slen=xmlStrlen(frag)+1;
      ax_sprintf(resource,0,slen,"%s",frag);
    }
  }
  else{
    frag=getProperty(node,(xmlChar*)"ID",(xmlChar*)"http://www.w3.org/1999/02/22-rdf-syntax-ns#");
    if(frag){
      slen=2+xmlStrlen(frag)+xmlStrlen(nodeBase);
      ax_sprintf(resource,0,slen,"%s#%s",nodeBase,frag);
    }
  }

  if(nodeBase)free(nodeBase);
  return retval;
}  

int RDFResource(xmlChar ***results, struct xpathObjWrap *objWrap){

  int i = 1, j = 0;
  xmlNodePtr node = NULL;
  xmlChar *resource = NULL;

  if ( !objWrap ) return(-1);
  if ( !results ) return(-1);

  locateNode(objWrap, i, &node, NULL, NULL);

  for (i = 2; node; i++){

      resourceURI(&resource,node);
      if(resource){
	*results=(xmlChar**)realloc(*results,(j+2)*sizeof(xmlChar*));
	*((*results)+j)=resource;
	resource=NULL;
	*((*results)+j+1)=NULL;
	j++;
      }

      locateNode(objWrap, i, &node, NULL, NULL);
  }
  return(0);
}


int RDFLiteral(xmlChar ***results, xmlChar *predicate, struct xpathObjWrap *objWrap){

  int i = 1, j = 0, slen = 0, retval = 0;
  xmlNodePtr node = NULL, cur = NULL;
  xmlChar* nodeURI = NULL;

  if (!objWrap) return(-1);
  if (!predicate) return(-1);
  if (!results) return(-1);

  nodeURI = (xmlChar*)malloc(1000 * sizeof(xmlChar));

  locateNode(objWrap, i, &node, NULL, NULL);

  for(i = 2; node; i++){
    
    node = node->children;
    nextSiblingOfType(&node,1);

    while(node){
    
      slen = xmlStrlen(node->name) + xmlStrlen(node->ns->href) + 1;
      ax_sprintf(&nodeURI, 1000, slen, "%s%s", node->ns->href, node->name);
    
      if (!xmlStrcmp(nodeURI, predicate)){
	    *results = (xmlChar**)realloc(*results, (j + 2) * sizeof(xmlChar*));
	    *((*results) + j + 1) = NULL;
	if ((cur = node->children)){
	  if (cur->type == XML_TEXT_NODE){
	    *((*results) + j) = xmlStrdup(cur->content);
	    j++;
	  }
	}
	else{
	  *((*results) + j) = (xmlChar*)malloc(sizeof(xmlChar));
	  **((*results) + j) = '\0';
	  j++;
	}
      }

      node = node->next;
      nextSiblingOfType(&node,1);
    }
    
    locateNode(objWrap, i, &node, NULL, NULL);
  }

  if (nodeURI) free(nodeURI);

  return retval;
}

int endResource(struct xpathObjWrap **objWrap, xmlChar *resource){

  int i = 1, k = 0;
  xmlChar *nodeResource = NULL;
  xmlNodeSetPtr matchingNodeSet = NULL;
  xmlXPathObjectPtr matchingObject = NULL;
  xmlNodePtr node = NULL;

  if (!resource) return(-1);
  if (!objWrap) return(-1);
  if (!*objWrap) return(0);

  locateNode((*objWrap), i, &node, NULL, NULL);
      
  for (i=2; node; i++){
    
    resourceURI(&nodeResource,node);
    
    if(nodeResource){
      if(!xmlStrcmp(nodeResource,resource)){
	
	if(!matchingNodeSet){
	  matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	  matchingNodeSet->nodeTab=NULL;
	}
	
	matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	matchingNodeSet->nodeTab[k]=node;
	k++;
      }
    }
    locateNode((*objWrap), i, &node, NULL, NULL);
  }
  
  if(matchingNodeSet){
    matchingObject=(xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));
    matchingObject->type=1;
    matchingObject->boolval=0;
    matchingObject->nodesetval=matchingNodeSet;
    matchingObject->nodesetval->nodeNr=k;
    matchingObject->nodesetval->nodeMax=k;
  }

  freeWrap(objWrap);
  *objWrap=NULL;
  wrapObject(objWrap,matchingObject);
  
  if(nodeResource)free(nodeResource);

  return(0);
}

int RDFTrace(struct xpathObjWrap **objWrap, xmlChar *predicate, int fb){

  /* fb == 0 -> trace forwards
     fb == 1 -> trace backwards
  */

  struct xpathObjWrap *result = NULL, *eval = NULL, *mWrap = NULL;
  int i = 1, j = 1, k = 0, retval = 0, slen = 0 ;
  xmlNodeSetPtr matchingNodeSet = NULL;
  xmlChar *fromResource = NULL, *toResource = NULL, *xptrExpr = NULL, *nodeBase = NULL;
  xmlChar *frag = NULL, *nodeURI = NULL, *rbase = NULL;
  xmlXPathObjectPtr matchingObject = NULL;
  xmlNodePtr node = NULL, cur = NULL;

  if (!objWrap) return(-1);
  if (!*objWrap) return(-1);
  if (!predicate) return(-1);
  if ( (fb != 0) && (fb != 1) ) return(-1);

  nodeURI = (xmlChar*)malloc(1000 * sizeof(xmlChar));
  xptrExpr = (xmlChar*)malloc(1000 * sizeof(xmlChar));
  toResource = (xmlChar*)malloc(1000 * sizeof(xmlChar));

  locateNode((*objWrap), i, &node, NULL, NULL);

  for(i = 2; node; i++){

    /* trace forwards */

    if(fb==0){

      slen=xmlStrlen(node->name)+xmlStrlen(node->ns->href)+1;
      ax_sprintf(&nodeURI,1000,slen,"%s%s",node->ns->href,node->name);
      
      if(!xmlStrcmp(predicate,(xmlChar*)"http://www.w3.org/2000/01/rdf-schema#type")){
	
	if(xmlStrcmp(nodeURI,(xmlChar*)"http://www.w3.org/1999/02/22-rdf-syntax-ns#Description")){

	  slen=128+(3*xmlStrlen(nodeURI));
	  ax_sprintf(&xptrExpr,1000,slen,"xpointer(//*[@rdf:ID=substring-after('%s',concat(xmlbase(),'#')) or @rdf:about=substring-after('%s',xmlbase()) or @rdf:about='%s'])",nodeURI,nodeURI,nodeURI);
	  
	  if(type_query(5,NULL,xptrExpr)) ERROR;
	  collectResults(&result);
	  if(type_query(6,NULL,xptrExpr)) ERROR;
	  collectResults(&result);
	  if(type_query(3,NULL,xptrExpr)) ERROR;
	  collectResults(&result);
	  
	}
      }
      
      node=node->children;
      nextSiblingOfType(&node,1);

      while(node){
	
	slen=xmlStrlen(node->name)+xmlStrlen(node->ns->href)+1;
	ax_sprintf(&nodeURI,1000,slen,"%s%s",node->ns->href,node->name);
	
	if(!xmlStrcmp(nodeURI,predicate)){
	  
	  frag=getProperty(node,(xmlChar*)"resource",(xmlChar*)"http://www.w3.org/1999/02/22-rdf-syntax-ns#");
	  
	  if(frag){

	    if(*frag=='#'){
	      rbase=xmlNodeGetBase(node->doc,node);
	      slen=xmlStrlen(rbase)+xmlStrlen(frag)+1;
	      ax_sprintf(&toResource,1000,slen,"%s%s",rbase,frag);
	      if (rbase) free (rbase);
	      rbase = NULL;
	    }

	    else{
	      slen=xmlStrlen(frag)+1;
	      ax_sprintf(&toResource,1000,slen,"%s",frag);
	    }
	    
	    slen=128+(3*xmlStrlen(toResource));
	    ax_sprintf(&xptrExpr,1000,slen,"xpointer(//*[@rdf:ID=substring-after('%s',concat(xmlbase(),'#')) or @rdf:about=substring-after('%s',xmlbase()) or @rdf:about='%s'])",toResource,toResource,toResource);
	    
	    if(type_query(5,NULL,xptrExpr)) ERROR;
	    collectResults(&result);
	    if(type_query(6,NULL,xptrExpr)) ERROR;
	    collectResults(&result);
	    if(type_query(3,NULL,xptrExpr)) ERROR;
	    collectResults(&result);
	    
	  }
	  
	  cur=node->children;
	  nextSiblingOfType(&cur,1);
	  
	  while(cur){
	    if(!matchingNodeSet){
	      matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	      matchingNodeSet->nodeTab=NULL;
	    }

	    matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	    matchingNodeSet->nodeTab[k]=cur;
	    k++;
	    cur=cur->next;
	    nextSiblingOfType(&cur,1);
	  }
	  /* add node itself, in case blank node */
	  if(!matchingNodeSet){
	    matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	    matchingNodeSet->nodeTab=NULL;
	  }
	  matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	  matchingNodeSet->nodeTab[k]=node;
	  k++;
	}
	node=node->next;
	nextSiblingOfType(&node,1);
      }
      
    }

    /* trace backwards */

    if(fb==1){

      resourceURI(&fromResource,node);

      if(!xmlStrcmp(predicate,(xmlChar*)"http://www.w3.org/2000/01/rdf-schema#type")){
	
	slen=55+xmlStrlen(fromResource);
	ax_sprintf(&xptrExpr,1000,slen,"xpointer(//*[concat(namespace-uri(),local-name())='%s'])",fromResource);
	
	if(type_query(5,NULL,xptrExpr)) ERROR;
	collectResults(&result);
	if(type_query(6,NULL,xptrExpr)) ERROR;
	collectResults(&result);
	if(type_query(3,NULL,xptrExpr)) ERROR;
	collectResults(&result);
      }
      
      if(fromResource){
	slen=80+xmlStrlen(fromResource);
	ax_sprintf(&xptrExpr,1000,slen,"xpointer(//*[@rdf:resource='%s' or @rdf:resource=substring-after('%s',xmlbase())])",fromResource,fromResource);
	
	if(eval){
	  freeWrap(&eval);
	  eval=NULL;
	}
	if(type_query(5,NULL,xptrExpr)) ERROR;
	collectResults(&eval);
	if(type_query(6,NULL,xptrExpr)) ERROR;
	collectResults(&eval);
	if(type_query(3,NULL,xptrExpr)) ERROR;
	collectResults(&eval);
	
	if(eval){
	  
	  locateNode(eval, j, &cur, NULL, NULL);
	    
	  for(j = 2; cur; j++){
	    
	    slen=xmlStrlen(cur->name)+xmlStrlen(cur->ns->href)+1;
	    ax_sprintf(&nodeURI,1000,slen,"%s%s",cur->ns->href,cur->name);
	    
	    if(!xmlStrcmp(nodeURI,predicate)){
	      if(!matchingNodeSet){
		matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
		matchingNodeSet->nodeTab=NULL;
	      }
	      matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	      matchingNodeSet->nodeTab[k]=cur->parent;
	      k++;
	    }
	    locateNode(eval, j, &cur, NULL, NULL);
	  }
	}
      }

      cur=node->parent;
      
      slen=xmlStrlen(cur->name)+xmlStrlen(cur->ns->href)+1;
      ax_sprintf(&nodeURI,1000,slen,"%s%s",cur->ns->href,cur->name);
      
      if(!xmlStrcmp(nodeURI,predicate)){
	if(!matchingNodeSet){
	  matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	  matchingNodeSet->nodeTab=NULL;
	}
	matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	matchingNodeSet->nodeTab[k]=cur->parent;
	k++;
      }

      /* add node's parent incase blank node */
      
      slen=xmlStrlen(node->name)+xmlStrlen(node->ns->href)+1;
      ax_sprintf(&nodeURI,1000,slen,"%s%s",node->ns->href,node->name);
      
      if(!xmlStrcmp(nodeURI,predicate)){
	if(!matchingNodeSet){
	  matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	  matchingNodeSet->nodeTab=NULL;
	}
	matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	matchingNodeSet->nodeTab[k]=node->parent;
	k++;
      }    
    }
    locateNode((*objWrap), i, &node, NULL, NULL);
  }
  
  if(matchingNodeSet){
    matchingObject=(xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));
    matchingObject->type=1;
    matchingObject->boolval=0;
    matchingObject->nodesetval=matchingNodeSet;
    matchingObject->nodesetval->nodeNr=k;
    matchingObject->nodesetval->nodeMax=k;
  }

  wrapObject(&mWrap,matchingObject);
  linkWrap(&result,mWrap);
  freeWrap(objWrap);
  *objWrap = result;

 cleanup:
  
  if (eval) freeWrap(&eval);
  if (xptrExpr) free(xptrExpr);
  if (fromResource) free(fromResource);
  if (toResource) free(toResource);
  if (nodeURI) free(nodeURI);
  if (nodeBase) free(nodeBase);
 
  return (0);
}

int findType(){

  int i = 1, retval = 0, *type_found = NULL, expHash = 0;
  xmlNodePtr node = NULL;
  xmlChar *resource = NULL;
  struct xpathObjWrap *objWrap = NULL;
  const xmlChar *name = NULL;
  struct cacheBtree* treeLocalPrev = NULL;
  int found = 0;
  
  if (gdata) resource = gdata->next_resource;
  else resource = ctxStrPtr->next_resource;

  expHash = toHash(resource);

  if (useCache == 1){
    
    searchCache(typeCacheOrigin, &treeLocalPrev, resource, (xmlChar*) "", (void**) &type_found, &found);

  }

  if (! found) {
    
    type_found = (int*) malloc (sizeof (int));
    *type_found = -1;
    
    startResource(&objWrap, resource);
    
    if (objWrap){
      
      locateNode(objWrap, i, &node, (NULL), (NULL));
      
      for (i = 2; node; i++){
	
	name = node->name;

	if (!xmlStrcmp(name, (xmlChar*)"Class")){
	  *type_found = 1;
	  break;
	}
	if (!xmlStrcmp(name, (xmlChar*)"DatatypeProperty")){
	  *type_found = 2;
	  break;
	}
	if (!xmlStrcmp(name,(xmlChar*)"ObjectProperty")){
	  *type_found = 3;
	  break;
	}
	if (!xmlStrcmp(name,(xmlChar*)"FunctionalProperty")){
	  *type_found = 4;
	  break;
	}
	
	locateNode(objWrap, i, &node, (NULL), (NULL));
	
      }
    }
    
    if (useCache == 1){
      
      addToCache(&treeLocalPrev, resource, (xmlChar*) "", (void*) type_found);
      if (!typeCacheOrigin) typeCacheOrigin = treeLocalPrev;
      
    }
  }
  
  if (objWrap) freeWrap(&objWrap);

  return *type_found;
  
}

int findCard(int* card, int* minCard, int* maxCard){

  xmlChar **ccard = NULL, **cminCard = NULL, **cmaxCard = NULL;
  int retval = 0;
  xmlChar *resource = NULL;
  struct xpathObjWrap *classObjects = NULL, *objWrap = NULL;

  if (gdata) resource = gdata->next_resource;
  else return (-1);

  startResource(&classObjects, gdata->class);
  linkWrap(&objWrap, classObjects);
  RDFTrace(&objWrap, (xmlChar*)"http://www.w3.org/2000/01/rdf-schema#range",0);
  linkWrap(&objWrap, classObjects);

  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2000/01/rdf-schema#subClassOf",0);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2002/07/owl#onProperty",0);
  endResource(&objWrap,resource);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2002/07/owl#onProperty",1);
  RDFLiteral(&ccard,(xmlChar*)"http://www.w3.org/2002/07/owl#cardinality",objWrap);
  RDFLiteral(&cminCard,(xmlChar*)"http://www.w3.org/2002/07/owl#minCardinality",objWrap);
  RDFLiteral(&cmaxCard,(xmlChar*)"http://www.w3.org/2002/07/owl#maxCardinality",objWrap);

  freeWrap(&objWrap);

  if (ccard) *card = atoi((const char*)*ccard);
  if (cminCard) *minCard = atoi ((const char*)*cminCard);
  if (cmaxCard) *maxCard = atoi ((const char*)*cmaxCard);

  return retval;
}


int isDelimiter(xmlChar *start, xmlChar *str, xmlChar **delimiter, int fb, int *count){

  int i=0,n=0;
  xmlChar* val=NULL;

  if(!count)return(-1);
  if(!str)return(-1);
  if(!start)return(-1);
  if(!delimiter)return(-1);
  if(str < start)return(-1);

  *count=0;

  for(i = 0;(val = delimiter[i]);i++){
    
    n = xmlStrlen(val); 
    
    if(fb == 1){
      
      if (*str == *val){
	
	if (!xmlStrncmp(str, val, n)){	    
	  *count = n;
	  return(0);
	}
      }
    }
    
    else {
      
      if((str - n + 1) >= start){
	
	if(*(str - n + 1) == *val){
	
	  if (!xmlStrncmp((str - n + 1),val,n)){
	    *count = - n;
	    return(0);
	  }
	}
      }
    }
  }

  return(0);
}

/* 
   Data element search function.
   Search forward (no > 0) or backward (no < 0) for the no(th) occurance
   of a string that IS NOT 'val' in the string 'orig'.
   The search starts at 'str'.
   Handles multiple delimiters.
   The start (begin) and end (end) of the string is marked.
*/

int ax_strstr(xmlChar *start, xmlChar *str, xmlChar **delimiter, int *sc, xmlChar** begin, xmlChar** end, int soff, int eoff) {

    int inc=1,counter=0,scinc=0,no=0,retval=0;
    xmlChar* strprev=NULL;
    xmlChar* smarker=NULL;
    xmlChar* emarker=NULL;
    
    if (!start) return(-1);
    if (!str) return(-1);
    if (!delimiter) return(-1);
    if (!sc) return(-1);

    no=*sc;

    if(*sc<0){
      inc=-1;
    }

    scinc=inc;

    if (begin) *begin = NULL;
    if (end) *end = NULL;

    while(str >= start && *str != '\0'){
      
      if(isDelimiter(start, str, delimiter, inc, &counter) == -1) return(-1);
      
      if(counter){
	
	if( !*sc && no ){
	  break;
	}
	
	while (counter) {
	  str+=inc;
	  counter-=inc;
	  if(*str == '\0' || str < start) break;
	}
	scinc=inc;

      }
      
      else {
	*sc-=scinc;
	if(scinc)strprev=str;
	str+=inc;
	scinc=0;
      }
    }
    
    if(!*sc){
      if(inc==1){
	if (end) *end = str - inc;
	if (begin) *begin = strprev;
      }
      else{
	if (begin) *begin = str - inc;
	if (end) *end = strprev;
      }
    }

    if (end && begin){
      if (soff < 0) smarker = *end + soff + 1;
      else if (soff > 0) smarker = *begin + soff - 1;
      if (smarker > *end || smarker < *begin) ERROR;
      
      if (eoff < 0) emarker = *end + eoff + 1;
      else if (eoff > 0) emarker = *begin + eoff - 1;
      if (emarker > *end || emarker < *begin) ERROR;
      
      *begin = smarker;
      *end = emarker;
    }

 cleanup:

    return retval;
}

int initItems(){

  int totalCount = 0, i = 0, retval = 0;
  int sindex = 0, eindex = 0, counter = 0;
  int inc = 1, imax = 0, icount = 0, noNodes = 0;
  xmlChar *resVars = NULL, *resExpr = NULL;
  xmlChar *begin = NULL, *end = NULL, **istr = NULL;
  xmlChar *cptr = NULL;

  if (!gdata) ERROR;
  if (!gdata->cdata_list) ERROR;
  if (!*(gdata->cdata_list)) ERROR;

  /* resolve startIndex and endIndex */

  if(gdata->startIndex){

    resolveVars(ctxStrPtr, *gdata->startIndex, &resVars);
    resolveExpr(resVars, &resExpr);
    sindex = atoi((const char*)resExpr);

  }

  if (sindex < 1) sindex = 1;

  if (gdata->endIndex){

    resolveVars(ctxStrPtr,*gdata->endIndex,&resVars);
    resolveExpr(resVars,&resExpr);
    eindex = atoi((const char*)resExpr);

  }

  if(eindex < sindex)eindex = 0;

  if(! gdata->delimiter){
    
    if( (noNodes = countNodes(gdata->item_set)) < sindex ) ERROR;
    
    gdata->charPtr = *( gdata->cdata_list + sindex - 1 );
    if (charData) { free (*charData); *charData = NULL; charData = NULL; }
    gdata->charData = gdata->charPtr;
    gdata->bufferPtr = gdata->charData;
    retval = noNodes;
    gdata->item_set_position = sindex;
    gdata->item_count = retval;
    gdata->index = sindex;
  }

  else{

    counter = (((sindex - gdata->index) - gdata->item_pos) + 1);
    if(counter < 0){
      inc = -1;
      counter--;
    }
    else counter++;

    istr = gdata->cdata_list;
    imax = countNodes(gdata->item_set);
    cptr = gdata->charPtr;

    for(i = gdata->item_set_position - 1; i > -1 && i < imax && counter ; i += inc){
      
      ax_strstr(istr[i],
		cptr,
		gdata->delimiter,
		&counter,
		&begin,
		&end,
		gdata->startChar,
		gdata->endChar);

      if ((i+inc) >  -1 && (i+inc) < imax) cptr = istr[i+inc];

    }

    gdata->item_set_position = ((i - inc) + 1);

    if (counter) ERROR;

    if(*(end + 1)){
      
      if ( charData ) { free (*charData); *charData = NULL; }
      gdata->charData = xmlStrndup(begin, ((end - begin) + 1));
      charData = &(gdata->charData);

      gdata->charPtr = begin;
      gdata->bufferPtr = gdata->charData;
      
    }
    
    else{
      
      gdata->charPtr = begin;
      if (charData) { free (*charData); *charData = NULL; charData = NULL; }
      gdata->charData = begin;
      gdata->bufferPtr = begin;
      
    }

    i-=inc;

    if (eindex) {

      icount = eindex - sindex + 1;
      counter = icount;
      cptr = gdata->charPtr;

      for(; i < imax && counter ; i ++){
	
	ax_strstr(istr[i],
		  cptr,
		  gdata->delimiter,
		  &counter,
		  (NULL),
		  (NULL),
		  gdata->startChar,
		  gdata->endChar);

	if ((i+1) < imax) cptr = istr[i+1];	
      }

      retval = icount - counter;
    }

      /* only count the remaining number of data elements when 
	 pushed into a corner!
      */      

    else{

      cptr = gdata->charPtr;

      for(; i < imax ; i ++){

	counter = 0;

	ax_strstr(istr[i],
		  cptr,
		  gdata->delimiter,
		  &counter,
		  (NULL),
		  (NULL),
		  gdata->startChar,
		  gdata->endChar);
	
	if ((i+1) < imax) cptr = istr[i+1];
	totalCount-=counter;

      }
      
      retval = totalCount;
      
    }
    
    gdata->item_count = retval;
    gdata->index = sindex;

  }

  if (gdata->property) free(gdata->property);
  gdata->property=gdata->next_resource;
  gdata->next_resource = NULL;

 cleanup:
  
  if(resVars)free(resVars);
  if(resExpr)free(resExpr);
  
  return(retval);
}

int Selectdi(){

  if (!gdata->item_set) return(-1);
  if (!gdata->item_set_position) return(-1);
  if (!gdata->charPtr) return(-1);
  if (!gdata->cdata_list) return(-1);

  gdata->set = &(gdata->item_set);
  gdata->set_pos = &(gdata->item_set_position);
  
  locateNode(gdata->item_set, gdata->item_set_position, &gdata->sel, (NULL), (NULL));
  gdata->sel_index = (int)((gdata->charPtr) - *(gdata->cdata_list + gdata->item_set_position - 1));

  if(gdata->charData){
    gdata->sel_index2 = gdata->sel_index + xmlStrlen(gdata->charData);
  }

  else {
    gdata->sel_index2 = gdata->sel_index + xmlStrlen(gdata->charPtr);
  }

  return (0);
}

int clearItems(){

  int i = 0;

  /* clear current cdata */
  
  if(gdata->cdata_list){

    if(gdata->freeable){

      for(i = 0; gdata->freeable[i]; i++){
	free(gdata->freeable[i]);
      }

    }

    free(gdata->cdata_list);
    gdata->cdata_list = NULL;
    gdata->freeable = NULL;
  }

  return(0);
}

int defaultXlocator(xmlChar ***xlocator, xmlChar *property, xmlChar *namespaceURI){

  int retval=0;

  if (!xlocator) ERROR;

  if((!*xlocator) && (supdocs == 0)){

    *xlocator = (xmlChar**)malloc(2*sizeof(xmlChar*));
    **(xlocator + 1) = NULL;
    **xlocator = (xmlChar*)malloc((127 + 3 * xmlStrlen(property) + xmlStrlen(namespaceURI)) * sizeof(xmlChar));
    sprintf((char*)(**xlocator),"xpointer(//*[local-name()='%s' or namespace-uri()='%s']|//eccp:complex[@role='http://www.grids.ac.uk/eccp/owl-ontologies#%s']|/*/@%s)",(char*)property,(char*)namespaceURI,(char*)property,(char*)property);

  }
  
 cleanup:

  return (retval);
}


int addNodeToObject(xmlXPathObjectPtr *object, xmlNodePtr node){

  int retval = 0, i = 0;
  xmlNodeSetPtr nodeSet = NULL;

  if (!object) ERROR
  if (!node) ERROR

  if (*object){
    nodeSet = (*object)->nodesetval;
  }
  
  else {

    *object = (xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));     
    (*object)->nodesetval = NULL;
    (*object)->type = 1;
    (*object)->boolval = 0;

  }

  if(!nodeSet){
    nodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
    nodeSet->nodeTab=NULL;
  }

  else i = nodeSet->nodeNr;

  nodeSet->nodeTab = realloc(nodeSet->nodeTab,(i + 2) * sizeof(xmlNodePtr));
  nodeSet->nodeTab[i] = node;
  i++;

  (*object)->nodesetval = nodeSet;
  (*object)->nodesetval->nodeNr = i;
  (*object)->nodesetval->nodeMax = i;

 cleanup:

  return (retval);
}

int checkPropertyRequest(){
  return (0);
}

int addToCache(struct cacheBtree** tree, xmlChar* expr1, xmlChar* expr2, void* load){

  struct cacheBtree* treeLocal = NULL;
  int expHash = 0, hash = 0;

  if (currentEvalCacheSize <= maxEvalCacheSize){

    expHash = toHash(expr1);
    
    treeLocal = (*tree);
    if (!(*tree)){
      treeLocal = (struct cacheBtree*)malloc(sizeof(struct cacheBtree));
      treeLocal->evalCachePtr = NULL;
      (*tree) = treeLocal;
      treeLocal->bigger = NULL;
      treeLocal->smaller = NULL;
    }
    else{
      hash = (*tree)->hash;
      if (expHash < hash){
	if (!(*tree)->smaller){
	  (*tree)->smaller = (struct cacheBtree*)malloc(sizeof(struct cacheBtree));
	  (*tree)->smaller->evalCachePtr = NULL;
	  (*tree)->smaller->smaller = NULL;
	  (*tree)->smaller->bigger = NULL;
	}
	treeLocal = (*tree)->smaller;
      }
      else if (expHash > hash){
	if (!(*tree)->bigger){
	  (*tree)->bigger = (struct cacheBtree*)malloc(sizeof(struct cacheBtree));
	  (*tree)->bigger->evalCachePtr = NULL;
	  (*tree)->bigger->smaller = NULL;
	  (*tree)->bigger->bigger = NULL;
	}
	treeLocal = (*tree)->bigger;
      }
    }
    
    if (!treeLocal->evalCachePtr){
      treeLocal->evalCacheOrigin = (struct evalCache*)malloc(sizeof(struct evalCache));
      treeLocal->evalCachePtr = treeLocal->evalCacheOrigin;
      treeLocal->evalCachePtr->prev = NULL;
    }
    else{
      treeLocal->evalCachePtr->next = (struct evalCache*)malloc(sizeof(struct evalCache));
      treeLocal->evalCachePtr->next->prev = treeLocal->evalCachePtr;
      treeLocal->evalCachePtr = treeLocal->evalCachePtr->next;
    }
    treeLocal->hash = expHash;
    treeLocal->evalCachePtr->next = NULL;
    treeLocal->evalCachePtr->expr2 = (xmlChar*)xmlStrdup(expr2);
    treeLocal->evalCachePtr->expr1 = (xmlChar*)xmlStrdup(expr1);
    treeLocal->evalCachePtr->load = (void*) load;
    currentEvalCacheSize ++;
  }
  return (0);
}


int searchCache(struct cacheBtree* tree, struct cacheBtree** finalBranch, xmlChar* expr1, xmlChar* expr2, void** load, int* found){

  struct evalCache* evalLocal = NULL;
  struct cacheBtree* branch = NULL;
  struct cacheBtree* branchPrev = NULL;
  int expHash = 0, hash = 0;

  if (!tree) return (-1);
  if (!expr1) return (-1);
  if (!expr2) return (-1);
  if (!load) return (-1);
  if (!finalBranch) return (-1);
  
  *load = NULL;
  branch = tree;
  branchPrev = branch;
  expHash = toHash(expr1);
  
  while(branch){
    hash = branch->hash;
    if(expHash < hash){
      branchPrev = branch;
      branch = branch->smaller;
    }
    else if(expHash > hash){
      branchPrev = branch;
      branch = branch->bigger;
    }
    else{
      evalLocal = branch->evalCacheOrigin;
      while (evalLocal){
	if (!xmlStrcmp(evalLocal->expr1, expr1)){
	  if (evalLocal->expr2){
	    if (!xmlStrcmp(evalLocal->expr2, expr2)){
	      *load = (void*) evalLocal->load;
	      *finalBranch = branchPrev;
	      *found = 1;
	      return (0);
	    }
	  }
	}
	evalLocal = evalLocal->next;
      }
      break;
    }
  }
  *finalBranch = branchPrev;
  return(0);
}

int getClassMap(struct classStr **classStrPtr, xmlChar *class){

  struct xpathObjWrap *objWrap = NULL;
  int expHash = 0, retval = 0 ;
  struct cacheBtree* treeLocalPrev = NULL;  
  xmlChar **dataSetType = NULL;
  int found = 0;
  
  expHash = toHash(class);

  if (useCache == 1){
    
    searchCache(classCacheOrigin, &treeLocalPrev, class, (xmlChar*) "", (void**) classStrPtr, &found);
    
  }
  
  if (! found){

    *classStrPtr = (struct classStr*) malloc (sizeof (struct classStr));
    (*classStrPtr)->xlocator = NULL;
    (*classStrPtr)->delimiter = NULL;
    (*classStrPtr)->dataSetType = 0;
    
    /* look for xlocator */
  
    startResource(&objWrap,class);
    RDFLiteral(&((*classStrPtr)->xlocator),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#xlocator",objWrap);

    /* look for dataSetType */
    
    RDFLiteral(&dataSetType,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#dataSetType",objWrap);

    if(dataSetType){
      if(xmlStrstr(*dataSetType,(xmlChar*)"virtual")){
	(*classStrPtr)->dataSetType = 1;

	/* look for delimiter */
	
	RDFLiteral(&((*classStrPtr)->delimiter),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#delimiter",objWrap);
	
      }
    }

    registerStandardDelimiters(&((*classStrPtr)->delimiter));
    
    if (useCache == 1){
      
      addToCache (&treeLocalPrev, class, (xmlChar*) "", (void*) *classStrPtr);
      if (!classCacheOrigin) classCacheOrigin = treeLocalPrev;
      
    }
  }
  
 cleanup:

  freeStringArray(&dataSetType);
  freeWrap(&objWrap);

  return retval;
}

int getPropMap(struct propertyStr **propertyStrPtr, xmlChar *property, xmlChar* class){

  struct xpathObjWrap *objWrap = NULL;
  xmlChar** charIndex = NULL;
  int expHash = 0, retval = 0 ;
  struct cacheBtree* treeLocalPrev = NULL;  
  int found = 0;
  
  expHash = toHash(property);

  if (useCache == 1){
    
    searchCache(propertyCacheOrigin, &treeLocalPrev, property, class, (void**) propertyStrPtr, &found);
    
  }
  
  if (! found){

    *propertyStrPtr = (struct propertyStr*) malloc (sizeof (struct propertyStr));
    (*propertyStrPtr)->xlocator = NULL;
    (*propertyStrPtr)->delimiter = NULL;
    (*propertyStrPtr)->sd = NULL;
    (*propertyStrPtr)->startIndex = NULL;
    (*propertyStrPtr)->endIndex = NULL;
    
    /* check global context */
    
    startResource(&objWrap, property);
    
    /* xlocator */
    RDFLiteral(&((*propertyStrPtr)->xlocator),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#xlocator",objWrap);
    
    /* delimiter */
    RDFLiteral(&((*propertyStrPtr)->delimiter),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#delimiter",objWrap);
    
    /* consecutive delimiter */
    RDFLiteral(&((*propertyStrPtr)->sd),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#consecutiveDelimiter",objWrap);
    
    /* startIndex */
    RDFLiteral(&((*propertyStrPtr)->startIndex),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#startIndex",objWrap);
    
    /* endIndex */
    RDFLiteral(&((*propertyStrPtr)->endIndex),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#endIndex",objWrap);
    
    /* startChar */
    RDFLiteral(&charIndex,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#startChar",objWrap);
    if(charIndex){
      (*propertyStrPtr)->startChar = atoi((char*)*charIndex);
      freeStringArray(&charIndex);
      charIndex = NULL;
    }
    else { (*propertyStrPtr)->startChar = 1; }
    
    /* endChar */
    RDFLiteral(&charIndex,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#endChar",objWrap);
    if(charIndex){
      (*propertyStrPtr)->endChar = atoi((char*)*charIndex);
      freeStringArray(&charIndex);
      charIndex = NULL;
    }
    else{ (*propertyStrPtr)->endChar = -1; }
    
    /* check class contexts */
    RDFTrace(&objWrap,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#context",0);
    RDFTrace(&objWrap,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#class",0);
    endResource(&objWrap,class);
    RDFTrace(&objWrap,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#class",1);
    
    /* xlocator */
    RDFLiteral(&((*propertyStrPtr)->xlocator),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#xlocator",objWrap);
    
    if (! (*propertyStrPtr)->xlocator) defaultXlocator(&((*propertyStrPtr)->xlocator), property, namespaceURI);
    
    if (! (*propertyStrPtr)->xlocator) END;
    
    /* delimiter */
    RDFLiteral(&((*propertyStrPtr)->delimiter),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#delimiter",objWrap);
    
    /* consecutive delimiter */
    RDFLiteral(&((*propertyStrPtr)->sd),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#consecutiveDelimiter",objWrap);
    
    /* startIndex */
    RDFLiteral(&((*propertyStrPtr)->startIndex),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#startIndex",objWrap);
    
    /* endIndex */
    RDFLiteral(&((*propertyStrPtr)->endIndex),(xmlChar*)"http://www.grids.ac.uk/eccp/ns#endIndex",objWrap);
    
    /* startChar */
    RDFLiteral(&charIndex,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#startChar",objWrap);
    if (charIndex){
      (*propertyStrPtr)->startChar = atoi((char*)*charIndex);
      free(charIndex);
      charIndex=NULL;
    }
    
    /* endChar */
    RDFLiteral(&charIndex,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#endChar",objWrap);
    if (charIndex) {
      (*propertyStrPtr)->endChar = atoi((char*)*charIndex);
      free(charIndex);
      charIndex=NULL;
    }
    
    registerStandardDelimiters(&((*propertyStrPtr)->delimiter));
    
    if (useCache == 1){
      
      addToCache (&treeLocalPrev, property, class, (void*) *propertyStrPtr);
      if (!propertyCacheOrigin) propertyCacheOrigin = treeLocalPrev;
      
    }
  }
  
 cleanup:
  
  freeWrap(&objWrap);
  return retval;
}

int copyStringArray(xmlChar*** newArrayPtr, xmlChar** arrayPtr) {
  return 0;;
}

int freeStringArray(xmlChar*** arrayPtr){

  int i = 0;
  
  if( *arrayPtr ){

      for(i = 0; *((*arrayPtr)+i) ; i++){
	free( *((*arrayPtr)+i) );
      }

      free ( *arrayPtr );
  }
  return (0);
}

int xml_get(xmlChar* property){

  xmlChar **sd = NULL;
  int i = 1;
  xmlChar **xlocator = NULL, **delimiter = NULL;
  xmlNodePtr context = NULL;
  xmlChar *resource = NULL;
  int retval = 0 ;
  struct xpathObjWrap *objWrap = NULL;
  struct xpathObjWrap *results = NULL;
  xmlXPathObjectPtr ro = NULL;
  struct propertyStr *propertyStrPtr = NULL;

  if (!ctxStrPtr) return (-1);
  if (!gdata) return (-1);
  if (!*(gdata->set)) return (-1);

  if (gdata->next_resource) free(gdata->next_resource);
  gdata->next_resource = (xmlChar*)malloc(((xmlStrlen(base_uri_string) + xmlStrlen(property)) + 1)*sizeof(xmlChar));
  resource = gdata->next_resource;
  sprintf((char*)resource,"%s%s",(char*)base_uri_string,(char*)property);

  /* check the request */

  if ( checkPropertyRequest() ) END;

  /* get property information */
  
  getPropMap ( &propertyStrPtr, resource, gdata->class);

  if ( propertyStrPtr ){
    
    xlocator = propertyStrPtr->xlocator;
    delimiter = propertyStrPtr->delimiter;
    sd = propertyStrPtr->sd;
    gdata->startChar = propertyStrPtr->startChar;
    gdata->endChar = propertyStrPtr->endChar;
    gdata->startIndex = propertyStrPtr->startIndex;
    gdata->endIndex = propertyStrPtr->endIndex;
    
  }

  /* all properties relating to virtual data sets must have a start index */

  /* set skip_delimiter */

  if(sd){
    if(xmlStrstr(*sd,(xmlChar*)"ignore")){
      gdata->skip_delimiter = 1;
    }
  }

  if ( ! xlocator ) ERROR;

  if ( **xlocator ) {

    /* get item set */
    
    /* if the physical data set is virtual, consider all the data sets */
    
    if (gdata->data_sets_type == 1){

      locateNode(*(gdata->set), i, &context, (NULL), (NULL));

      for (i = 2; context; i++){

	type_query(5,context,*xlocator);
	collectResults(&results);
	locateNode(*(gdata->set), i, &context, (NULL), (NULL));	
	
      }
    }
    
    else{
      
      context = gdata->sel;
      type_query(5, context, *xlocator);
      collectResults(&results);

      gdata->delimiter = delimiter;
      
    }

    if ( results ){
      if ( gdata->item_set ) freeWrap ( &gdata->item_set );
      gdata->item_set = results;
    }
    else END
    
    clearItems();
    getDataElements( gdata->item_set, &(gdata->cdata_list), (NULL), (NULL));
    gdata->charPtr = gdata->cdata_list[0];
    gdata->item_set_position = 1;
    gdata->item_pos = 1;
    gdata->index = 1;
    
  }
  
  else{
    
    if ( gdata->data_sets_type == 1 ){
      if ( !gdata->item_set ){

	gdata->item_set = *(gdata->set);
	incWrapRefs (gdata->item_set, 1);

	clearItems();
	getDataElements( gdata->item_set, &(gdata->cdata_list), (NULL), (NULL) );
	gdata->charPtr = gdata->cdata_list[0];
	gdata->item_set_position = 1;
	gdata->index = 1;
      }
      gdata->item_pos = 1;
    }
    
    else{
      
      addNodeToObject(&ro,gdata->sel);
      wrapObject(&results,ro);
      gdata->item_set = results;

      clearItems();
      getDataElements(gdata->item_set, &(gdata->cdata_list), (NULL), &(gdata->freeable));

      gdata->delimiter = delimiter;
      gdata->charPtr = gdata->cdata_list[0];
      gdata->item_set_position = 1;
      gdata->item_pos = 1;
      gdata->index = 1;
      
    }
    
  }

  retval = initItems();
  if (retval > 0){
    Selectdi();
  }
  
 cleanup:
  
  return retval;
  
}

xmlChar* getProperty(xmlNodePtr node,xmlChar* property,xmlChar* ns){

  xmlAttrPtr propLocal=NULL;
  xmlChar* result=NULL;

  if(!node)return (NULL);
  if(!property)return (NULL);

  propLocal=node->properties;
  while(propLocal){
    if(propLocal->name){
      if(!xmlStrcmp(propLocal->name,property)){
	if(propLocal->ns){
	  if(propLocal->ns->href){
	    if(!xmlStrcmp(propLocal->ns->href,ns)){
	      result=propLocal->children->content;
	      break;
	    }
	  }
	}
	else if(!ns){
	  result=propLocal->children->content;
	  break;
	}
      }
    }
    propLocal=propLocal->next;
  }
  
  return result;
}

int getNodes(int inst, xmlChar* exp, struct ctxStr* ctxStrLocal, struct xpathObjWrap **nodes){

  xmlChar *resVars = NULL, *resExpr = NULL;
  int l = 0, retval = 0;
  struct ctxStr* ctxStrLocal2 = NULL;

  for ( l = 1; l <= inst; l++ ){

    sprintf((char*)ctxStrLocal->var2, "%d", l);
    resolveVars(ctxStrLocal, exp, &resVars);
    resolveExpr(resVars, &resExpr);
    ctxStrLocal2 = ctxStrOrigin;

    while (ctxStrLocal2){
      if
	
	(((ctxStrLocal->type == 1) 
	  && ((ctxStrLocal2->parent == ctxStrLocal) 
	      || (ctxStrLocal2 == ctxStrLocal))) 
	 || ((ctxStrLocal->type == 2) 
	     && ((ctxStrLocal2 == ctxStrLocal->parent) 
		 || (ctxStrLocal2->parent == ctxStrLocal->parent))) 
	 || ((ctxStrLocal->type == 3) 
	     && (ctxStrLocal2->type != 3)))
	
	{
	
	if ( xptr_query(ctxStrLocal2, ctxStrLocal2->original_doc_root, resExpr) ) ERROR;
	linkWrap(nodes, ctxStrLocal2->xpathObj);
	freeWrap(&ctxStrLocal2->xpathObj);
	ctxStrLocal2->xpathObj = NULL;
      }
      ctxStrLocal2 = ctxStrLocal2->next;
    }
    if (resVars) free(resVars);
    resVars = NULL;
    if (resExpr) free(resExpr);
    resExpr = NULL;
  }

 cleanup:

  return retval;
}

int resolveLink(xmlNodePtr node, struct ctxStr* ctxStrLocal){

  int scope = 0;
  int j = 1, k = 1, l = 0, retval = 0, len = 0, linkinstances, locinstances;
  xmlNodePtr context = NULL, fromNode = NULL, arcNode = NULL;
  xmlNodePtr fromLocatorNode = NULL, toLocatorNode = NULL;
  struct xpathObjWrap *arcNodes = NULL, *fromNodes = NULL;
  struct xpathObjWrap *toNodes = NULL, *fromLocatorNodes = NULL, *toLocatorNodes =NULL;
  struct xpathObjWrap *commonToNodes = NULL, *commonFromNodes = NULL;
  xmlChar *arcFromVal = NULL, *arcToVal = NULL, *xptrExpr = NULL, *inst = NULL, *scopestr = NULL;
  xmlChar* idhrefValue = NULL;
  int start = -1, end = -1;

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));

  /* get value of the instances attribute */
  
  inst = getProperty(node,(xmlChar*)"instances",(NULL));
  
  scopestr = getProperty(node,(xmlChar*)"scope",(NULL));
  
  /* if instances attribute not found, default to 1 */
  
  if (!inst) linkinstances = 1;
  else {
    linkinstances = atoi((const char*)inst);
    inst = NULL;
  }
  
  /* if scope attribute not found, default to 0 */
  
  if (!scopestr) scope = 0;
  else {
    scope = atoi((const char*)scopestr);
    scopestr=NULL;
  }
  
  /* look for arc nodes */
  
  context = node;
  sprintf((char*)xptrExpr, "xpointer(//eccp:arc | //molpro:associate)");
  if (xptr_query(ctxStrPtr,context,xptrExpr)) ERROR;
  if (arcNodes) freeWrap (&arcNodes);
  linkWrap(&arcNodes, ctxStrPtr->xpathObj);
  freeWrap(&ctxStrPtr->xpathObj);
  ctxStrPtr->xpathObj=NULL;
  
  if(arcNodes){
    
    j=1;
    
    locateNode(arcNodes, j, &arcNode, (NULL), (NULL));
		
    for(j = 2; arcNode; j++){
      
      /* get the xlink:from and xlink:to attribute values */
      
      arcFromVal=getProperty(arcNode,(xmlChar*)"from",(xmlChar*)"http://www.w3.org/1999/xlink");
      arcToVal=getProperty(arcNode,(xmlChar*)"to",(xmlChar*)"http://www.w3.org/1999/xlink");
      
      /* get from locator nodes */
      
      if(arcFromVal){
	len=xmlStrlen(arcFromVal)+xmlStrlen((const xmlChar*)"xpointer(//eccp:locator[@xlink:label='']|//molpro:atoms[@xlink:label=''])");
	if(len>99){
	  xptrExpr=realloc(xptrExpr,(len+1)*sizeof(xmlChar));
	}
	sprintf((char*)xptrExpr,"xpointer(//eccp:locator[@xlink:label='%s']|//molpro:atoms[@xlink:label='%s'])",(char*)arcFromVal,(char*)arcFromVal);
      }
      
      else{
	sprintf((char*)xptrExpr,"xpointer(//eccp:locator|//molpro:atoms)");
      }
      
      context = node;
      if (xptr_query(ctxStrPtr, context, xptrExpr)) ERROR;
      if (fromLocatorNodes) freeWrap (&fromLocatorNodes);
      linkWrap(&fromLocatorNodes, ctxStrPtr->xpathObj);
      freeWrap(&ctxStrPtr->xpathObj);
      ctxStrPtr->xpathObj = NULL;
      
      /* get to locator nodes */
      
      if(arcToVal){
	len=xmlStrlen(arcToVal)+xmlStrlen((const xmlChar*)"xpointer(//eccp:locator[@xlink:label='']|//molpro:bases[@xlink:label=''])");
	if(len>99){
	  xptrExpr=realloc(xptrExpr,(len+1)*sizeof(xmlChar));
	}
	sprintf((char*)xptrExpr,"xpointer(//eccp:locator[@xlink:label='%s']|//molpro:bases[@xlink:label='%s'])",(char*)arcToVal,(char*)arcToVal);
      }
      else{
	sprintf((char*)xptrExpr,"xpointer(//eccp:locator|//molpro:bases)");
      }
      
      if(xptr_query(ctxStrPtr,context,xptrExpr)) ERROR;
      if (toLocatorNodes) freeWrap (&toLocatorNodes);
      linkWrap(&toLocatorNodes,ctxStrPtr->xpathObj);
      freeWrap(&ctxStrPtr->xpathObj);
      ctxStrPtr->xpathObj=NULL;
      
      /* complete common evaluations */
      
      k = 1;
      
      locateNode(fromLocatorNodes, k, &fromLocatorNode, (NULL), (NULL));
      
      for(k = 2; fromLocatorNode; k++){
	
	/* get the xlink:href and instances attribute values */ 
	
	idhrefValue=getProperty(fromLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");

	if (*idhrefValue == '#') idhrefValue ++;
	
	inst=getProperty(fromLocatorNode,(xmlChar*)"instances",(NULL));
	
	if(!xmlStrstr(idhrefValue,(const xmlChar*)"$1")){
	  
	  /* if instances attribute not found, default to 1 */
	  
	  if(!inst)locinstances=1;
	  else {
	    sscanf((char*)inst,"%d",&locinstances);
	    inst=NULL;
	  }
	  
	  if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(xmlChar));
	  
	  /* get the common from nodes */
	  /* queries are across all documents */
	  
	  getNodes(locinstances,idhrefValue,ctxStrLocal,&commonFromNodes);
	  
	  idhrefValue=NULL;
	}
	locateNode(fromLocatorNodes, k, &fromLocatorNode, (NULL), (NULL));
      }
      
      k = 1;
      
      locateNode(toLocatorNodes, k, &toLocatorNode, (NULL), (NULL));
      
      for(k = 2; toLocatorNode; k++){
	
	/* get the xlink:href and instances attribute values */
	
	idhrefValue=getProperty(toLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");

	if (*idhrefValue == '#') idhrefValue ++;
	
	inst=getProperty(toLocatorNode,(xmlChar*)"instances",(NULL));
	
	if(!xmlStrstr((const xmlChar*)idhrefValue,(const xmlChar*)"$1")){
	  
	  if(!inst)locinstances=1;
	  else {
	    sscanf((char*)inst,"%d",&locinstances);
	    inst=NULL;
	  }
	  
	  if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(char));
	  
	  /* get the common to nodes */
	  /* queries are across all documents */
	  
	  
	  getNodes(locinstances,idhrefValue,ctxStrLocal,&commonToNodes);
	  idhrefValue=NULL;
	}
	locateNode(toLocatorNodes, k, &toLocatorNode, (NULL), (NULL));
      }
      
      if (!ctxStrLocal->var) ctxStrLocal->var = (xmlChar*)malloc(100*sizeof(char));
      
      for( k = 1; k <= linkinstances; k++){
	
	sprintf((char*)ctxStrLocal->var, "%d", k);
	
	l = 1;
	
	locateNode(fromLocatorNodes, l, &fromLocatorNode, (NULL), (NULL));
	
	for(l = 2; fromLocatorNode; l++){
	  
	  /* get the xlink:href and instances attribute values */
	  
	  idhrefValue=getProperty(fromLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");
	  
	  inst=getProperty(fromLocatorNode,(xmlChar*)"instances",(NULL));
	  
	  /* if instances attribute not found, default to 1 */
	  
	  if(xmlStrstr((const xmlChar*)idhrefValue,(xmlChar*)"$1")){
	    
	    if(!inst)locinstances=1;
	    else {
	      sscanf((char*)inst,"%d",&locinstances);
	      inst=NULL;
	    }
	    
	    fromNodes=NULL;
	    if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(char));
	    
	    /* get the from nodes */
	    /* queries are across all documents */
	    
	    getNodes(locinstances,idhrefValue,ctxStrLocal,&fromNodes);
	    
	    idhrefValue=NULL;
	  }
	  locateNode(fromLocatorNodes, l, &fromLocatorNode, (NULL), (NULL));
	}
	
	l = 1;
	
	locateNode(toLocatorNodes, l, &toLocatorNode, (NULL), (NULL));
	
	for(l = 2; toLocatorNode; l++){
	  
	  /* get the xlink:href and instances attribute values */
	  
	  idhrefValue=getProperty(toLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");
	  
	  inst=getProperty(toLocatorNode,(xmlChar*)"instances",(NULL));
	  
	  if(xmlStrstr(idhrefValue,(xmlChar*)"$1")){
	    
	    if(!inst)locinstances=1;
	    else {
	      sscanf((char*)inst,"%d",&locinstances);
	      inst=NULL;
	    }
	    
	    toNodes=NULL;
	    if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(xmlChar));
	    
	    /* get the to nodes */
	    /* queries are across all documents */
	    
	    getNodes(locinstances,idhrefValue,ctxStrLocal,&toNodes);
	    
	    idhrefValue=NULL;
	  }
	  locateNode(toLocatorNodes, l, &toLocatorNode, (NULL), (NULL));
	}
	
	linkWrap(&fromNodes,commonFromNodes);
	linkWrap(&toNodes,commonToNodes);
	
	if((fromNodes)&&(toNodes)){
	  
	  l = 1;
	  locateNode(fromNodes, l, &fromNode, &start, &end);
	  
	  for(l = 2; fromNode; l++){
		      
	    if(fromNode->_private == NULL){
	      createprivate(fromNode);
	      ((rsPtr)fromNode->_private)->linksPtr=(lsPtr)malloc(sizeof(links));
	      ((rsPtr)fromNode->_private)->linksOrigin=((rsPtr)fromNode->_private)->linksPtr;
	    }
	    else{
	      ((rsPtr)fromNode->_private)->linksPtr->next=(lsPtr)malloc(sizeof(links));
	      ((rsPtr)fromNode->_private)->linksPtr->next->prev=((rsPtr)fromNode->_private)->linksPtr;
	      ((rsPtr)fromNode->_private)->linksPtr=((rsPtr)fromNode->_private)->linksPtr->next;
	    }
	    ((rsPtr)fromNode->_private)->linksPtr->next=NULL;
	    ((rsPtr)fromNode->_private)->linksPtr->nodes=NULL;
	    ((rsPtr)fromNode->_private)->linksPtr->start=start;
	    ((rsPtr)fromNode->_private)->linksPtr->end=end;
	    ((rsPtr)fromNode->_private)->linksPtr->scope=scope;
	    
	    linkWrap(&(((rsPtr)fromNode->_private)->linksPtr->nodes),toNodes);
	    ((rsPtr)fromNode->_private)->linksPtr->key=((rsPtr)node->_private)->copies[k-1];
	    locateNode(fromNodes, l, &fromNode, &start, &end);
	  }
	}
	
	if(toNodes){
	  freeWrap(&toNodes);
	  toNodes=NULL;
	}
	
	if(fromNodes){
	  freeWrap(&fromNodes);
	  fromNodes=NULL;
	}
      }
      
      if(commonToNodes){
	freeWrap(&commonToNodes);
	commonToNodes=NULL;
      }
      
      if(commonFromNodes){
	freeWrap(&commonFromNodes);
	commonFromNodes=NULL;
      }
      locateNode(arcNodes, j, &arcNode, (NULL), (NULL));
    }
  }

 cleanup:

  if (fromLocatorNodes) freeWrap(&fromLocatorNodes);
  if (toLocatorNodes) freeWrap(&toLocatorNodes);
  if (arcNodes) freeWrap(&arcNodes);
  if (xptrExpr) free (xptrExpr);

  return (0);
}

int resolveAssociate(xmlNodePtr node, struct ctxStr* ctxStrLocal){

  int j = 1, retval = 0;
  xmlNodePtr context = NULL, fromNode = NULL;
  struct xpathObjWrap *fromNodes = NULL, *toNodes = NULL;
  xmlChar *xptrExpr = NULL;

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));

  /* look for children of the from nodes */
  
  context=node;
  sprintf((char*)xptrExpr,"xpointer(//eccp:from/child::*)");
  if(xptr_query(ctxStrPtr,context,xptrExpr)) ERROR;
  linkWrap(&fromNodes,ctxStrPtr->xpathObj);
  freeWrap(&ctxStrPtr->xpathObj);
  ctxStrPtr->xpathObj=NULL;
  
  /* look for children of the to nodes */
  
  sprintf((char*)xptrExpr,"xpointer(//eccp:to/child::*)");
  if(xptr_query(ctxStrPtr,context,xptrExpr)) ERROR;
  linkWrap(&toNodes,ctxStrPtr->xpathObj);
  freeWrap(&ctxStrPtr->xpathObj);
  ctxStrPtr->xpathObj=NULL;
  
  if((fromNodes)&&(toNodes)){
    
    j = 1;
    locateNode(fromNodes,j,&fromNode,(NULL),(NULL));
    
    for(j = 2; fromNode; j++){
      
      if(fromNode->_private==NULL){
	createprivate(fromNode);
	((rsPtr)fromNode->_private)->linksPtr=(lsPtr)malloc(sizeof(links));
	((rsPtr)fromNode->_private)->linksOrigin=((rsPtr)fromNode->_private)->linksPtr;
      }
      else{
	((rsPtr)fromNode->_private)->linksPtr->next=(lsPtr)malloc(sizeof(links));
	((rsPtr)fromNode->_private)->linksPtr->next->prev=((rsPtr)fromNode->_private)->linksPtr;
	((rsPtr)fromNode->_private)->linksPtr=((rsPtr)fromNode->_private)->linksPtr->next;
      }
      ((rsPtr)fromNode->_private)->linksPtr->next=NULL;
      ((rsPtr)fromNode->_private)->linksPtr->nodes=NULL;
      ((rsPtr)fromNode->_private)->linksPtr->start=-1;
      ((rsPtr)fromNode->_private)->linksPtr->end=-1;
      
      linkWrap(&(((rsPtr)fromNode->_private)->linksPtr->nodes),toNodes);
      ((rsPtr)fromNode->_private)->linksPtr->key=node;
      locateNode(fromNodes, j, &fromNode, (NULL), (NULL));
    }
  }
  if(toNodes){
    freeWrap(&toNodes);
    toNodes=NULL;
  }
  
  if(fromNodes){
    freeWrap(&fromNodes);
    fromNodes=NULL;
  }

 cleanup:

  return (0);
}

int resolveComplex(xmlNodePtr node, struct ctxStr* ctxStrLocal){

  int j = 1, k = 1, l = 0, retval = 0, linkinstances, locinstances;
  xmlNodePtr context = NULL, fromNode = NULL, toNode = NULL;
  xmlNodePtr toLocatorNode = NULL;
  struct xpathObjWrap *fromNodes = NULL;
  struct xpathObjWrap *toNodes = NULL, *toLocatorNodes =NULL;
  struct xpathObjWrap *commonToNodes = NULL, *commonFromNodes = NULL;
  xmlChar *xptrExpr = NULL, *inst = NULL;
  xmlChar* idhrefValue = NULL;
  xmlXPathObjectPtr fromNodesObj = NULL;

  xptrExpr = (xmlChar*)malloc(100 * sizeof(xmlChar));

  /* get value of the instances attribute */
  
  inst = getProperty(node, (xmlChar*)"instances", (NULL));
  
  /* if instances attribute not found, default to 1 */
  
  if(!inst)linkinstances=1;
  else {
    sscanf((char*)inst,"%d",&linkinstances);
    inst=NULL;
  }
  
  /* find the locator nodes */
  
  context=node;
  sprintf((char*)xptrExpr,"xpointer(//eccp:locator)");
  if(xptr_query(ctxStrPtr,context,xptrExpr)) ERROR;
  linkWrap(&toLocatorNodes,ctxStrPtr->xpathObj);
  freeWrap(&ctxStrPtr->xpathObj);
  ctxStrPtr->xpathObj=NULL;
  
  /* complete common evaluations */
  
  if (toLocatorNodes){
    
    j = 1;
    
    locateNode(toLocatorNodes, j, &toLocatorNode, (NULL), (NULL));
    
    for( j = 2; toLocatorNode; j++){
      
      /* get the xlink:href and instances attribute values */
      
      idhrefValue=getProperty(toLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");
      
      inst=getProperty(toLocatorNode,(xmlChar*)"instances",(NULL));
      
      if(!xmlStrstr(idhrefValue,(xmlChar*)"$1")){
	
	/* if instances attribute not found, default to 1 */
	
	if(!inst)locinstances=1;
	else {
	  sscanf((char*)inst,"%d",&locinstances);
	  inst=NULL;
	}
	
	if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(xmlChar));
	
	getNodes(locinstances,idhrefValue,ctxStrLocal,&commonToNodes);
	
      }
      locateNode(toLocatorNodes, j, &toLocatorNode, (NULL), (NULL));
    }
  }
  
  /* complete non-common evaluations */
  
  if(!ctxStrLocal->var)ctxStrLocal->var=malloc(100*sizeof(xmlChar));
  
  for(j = 1; j <= linkinstances; j++){
    
    sprintf((char*)ctxStrLocal->var,"%d",j);
    
    if(toLocatorNodes){
      
      k = 1;
      locateNode(toLocatorNodes,k,&toLocatorNode,(NULL),(NULL));
      
      for(k = 2; toLocatorNode; k++){
	
	/* get the xlink:href and instances attribute values */
	
	idhrefValue=getProperty(toLocatorNode,(xmlChar*)"href",(xmlChar*)"http://www.w3.org/1999/xlink");
	
	inst=getProperty(toLocatorNode,(xmlChar*)"instances",(NULL));
	
	/* if instances attribute not found, default to 1 */
	
	if(xmlStrstr(idhrefValue,(xmlChar*)"$1")){
	  
	  if(!inst)locinstances=1;
	  else {
	    sscanf((char*)inst,"%d",&locinstances);
	    inst=NULL;
	  }
	  
	  if(!ctxStrLocal->var2)ctxStrLocal->var2=(xmlChar*)malloc(100*sizeof(xmlChar));
	  
	  getNodes(locinstances,idhrefValue,ctxStrLocal,&toNodes);
	}
	locateNode(toLocatorNodes, k, &toLocatorNode, (NULL), (NULL));
      }
      
      fromNodesObj=(xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));
      fromNodesObj->nodesetval=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
      fromNodesObj->nodesetval->nodeNr=1;
      fromNodesObj->boolval=0;
      fromNodesObj->type=1;
      fromNodesObj->nodesetval->nodeTab=(xmlNodePtr*)malloc(sizeof(xmlNodePtr));
      fromNodesObj->nodesetval->nodeTab[0]=((rsPtr)node->_private)->copies[j-1];
      
      wrapObject(&fromNodes,fromNodesObj);
      fromNode=fromNodesObj->nodesetval->nodeTab[0];
      
      linkWrap(&toNodes,commonToNodes);
      
      if((fromNodes)&&(toNodes)){
	
	l = 1;
	locateNode(toNodes, l, &toNode, (NULL), (NULL));
	
	for(l = 2; toNode; l++){
	  
	  if(toNode->_private==NULL){
	    createprivate(toNode);
	  }
	  linkWrap(&(((rsPtr)toNode->_private)->associations),fromNodes);
	  locateNode(toNodes, l, &toNode, (NULL), (NULL));
	}
	
	if(fromNode->_private==NULL){
	  createprivate(fromNode);				      
	}
	linkWrap(&(((rsPtr)fromNode->_private)->local),toNodes);
      }
    }
    
    if(toNodes){
      freeWrap(&toNodes);
      toNodes=NULL;
    }
    
    if(fromNodes){
      freeWrap(&fromNodes);
      fromNodes=NULL;
    }
    
    if(commonFromNodes){
      freeWrap(&commonFromNodes);
      commonFromNodes=NULL;
    }
  }
  if(commonToNodes){
    freeWrap(&commonToNodes);
    commonToNodes=NULL;
  }

 cleanup:
  
  if (toLocatorNodes) freeWrap(&toLocatorNodes);
  if (xptrExpr) free (xptrExpr);
  return(0);
}

int resolveLinks(){

  int i = 1, retval = 0;
  xmlNodePtr context = NULL, node = NULL;
  struct xpathObjWrap *nodes = NULL, *arcNodes = NULL;
  struct xpathObjWrap *fromLocatorNodes = NULL, *toLocatorNodes =NULL;
  xmlChar *xptrExpr = NULL;
  struct ctxStr* ctxStrLocal = NULL;

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));
  ctxStrLocal=ctxStrOrigin;

  ctxStrLocal=ctxStrOrigin;

  while(ctxStrLocal){
      
    /* find associate nodes */

    sprintf((char*)xptrExpr,"xpointer(//eccp:associate)");
    context=ctxStrLocal->original_doc_root;
    if(xptr_query(ctxStrLocal,context,xptrExpr)) ERROR;
    linkWrap(&nodes,ctxStrLocal->xpathObj);
    freeWrap(&ctxStrLocal->xpathObj);
    ctxStrLocal->xpathObj=NULL;

    /* add original complex and link nodes */

    linkWrap(&nodes,ctxStrLocal->copyNodes);

    i = 1;

    locateNode(nodes, i, &node, (NULL), (NULL));

    for (i = 2; node; i++){
      
      if ( !xmlStrcmp(node->ns->href, (xmlChar*)"http://www.grids.ac.uk/eccp/ns#") ){
	
	if ( !xmlStrcmp(node->name, (xmlChar*)"link") ) resolveLink(node,ctxStrLocal);	
	if ( !xmlStrcmp((xmlChar*)node->name, (xmlChar*)"associate") ) resolveAssociate(node,ctxStrLocal);
	if ( !xmlStrcmp(node->name, (xmlChar*)"complex") ) resolveComplex(node,ctxStrLocal);
	
      }

      if ( !xmlStrcmp(node->ns->href, (xmlChar*)"http://www.molpro.net/schema/molpro2005") ){
	
	if ( !xmlStrcmp(node->name, (xmlChar*)"association") ) resolveLink(node,ctxStrLocal);	
	
      }

      locateNode(nodes, i, &node, (NULL), (NULL)); 
    }

    freeWrap(&nodes);
    nodes = NULL;
    ctxStrLocal = ctxStrLocal->next;
  }
  
  linksResolved = 1;
  
 cleanup:
  
  if (arcNodes) freeWrap(&arcNodes);
  if (fromLocatorNodes) freeWrap(&fromLocatorNodes);
  if (toLocatorNodes) freeWrap(&toLocatorNodes);
  if (xptrExpr) free(xptrExpr);

  /* need to clear up common to and from nodes */

  return (retval);
  
}

void message_handler(xmlChar* function_name,xmlChar* message,int type){   /* error handling routine */
 
  /* to do */

  printf("%s returned: %s\n",(char*)function_name,(char*)message);
  type=0;
}

int getKeys(struct xpathObjWrap** keys, struct xpathObjWrap* objects){

  int i = 1;
  xmlNodePtr node=NULL;
  xmlXPathObjectPtr keysObj=NULL;

  if(!objects)return(-1);

  locateNode(objects, i, &node, (NULL), (NULL));
  
  for( i = 2; node; i++){

    if(((!xmlStrcmp(node->ns->href,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#"))&&(!xmlStrcmp(node->name,(xmlChar*)"link")))||((!xmlStrcmp(node->ns->href,(xmlChar*)"http://www.molpro.net/schema/molpro2005"))&&(!xmlStrcmp(node->name,(xmlChar*)"association")))){
      if(!keysObj){
	keysObj=(xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));
	keysObj->nodesetval=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	keysObj->nodesetval->nodeTab=NULL;
	keysObj->nodesetval->nodeNr=0;
      }
  
      keysObj->nodesetval->nodeTab=realloc(keysObj->nodesetval->nodeTab,(keysObj->nodesetval->nodeNr+1)*sizeof(xmlNodePtr));
      keysObj->nodesetval->nodeTab[keysObj->nodesetval->nodeNr]=node;
      keysObj->nodesetval->nodeNr++;
    }
    locateNode(objects, i, &node, (NULL), (NULL));
  }

  if(keysObj){

    keysObj->nodesetval->nodeTab=realloc(keysObj->nodesetval->nodeTab,(keysObj->nodesetval->nodeNr+1)*sizeof(xmlNodePtr));
    keysObj->nodesetval->nodeTab[keysObj->nodesetval->nodeNr]=NULL;
    
    keysObj->boolval=0;
    keysObj->type=1;
    
    wrapObject(keys,keysObj);
  }
  return (0);
}

int getevalcontext(xmlNodePtr context,xmlChar** evalContext){

  int i=0,size=0,newsize=0,ecsize=0,ectmpsize=0;
  xmlChar* ectmp=NULL;
  struct linkupinfo* linkinfo=NULL;
  xmlNodePtr node=NULL;

  if(!context)return(-1);
  if(!evalContext)return(-1);
  
  ectmp=(xmlChar*)malloc(100*sizeof(xmlChar));
  *evalContext=(xmlChar*)malloc(100*sizeof(xmlChar));
  size=100;
  sprintf((char*)*evalContext,"%d",(int)context);

  if(context->_private){
    linkinfo=((rsPtr)context->_private)->linkinfo;
  }

  if(!linkinfo){goto cleanup;}
  
  node=context->last;
  
  for(i=0;i<linkinfo->no;i++){
    
    sprintf((char*)ectmp,"%d",(int)node);
    ectmpsize=xmlStrlen(ectmp);
    ecsize=xmlStrlen(*evalContext);
    newsize=ectmpsize+ecsize;
    
    while(newsize>=size){
      *evalContext=realloc(*evalContext,(size+100)*sizeof(xmlChar));
      size+=100;
    }
    sprintf((char*)((*evalContext)+ecsize),"%d",(int)node);
    node=node->prev;
  }

 cleanup:
  
  if(ectmp)free(ectmp);

  return(0);
}


int createlinkinfo(xmlNodePtr currentNode){

  struct linkupinfo* linkinfo=NULL;

  if(!currentNode->_private)return(-1);

  linkinfo=(struct linkupinfo*)malloc(sizeof(struct linkupinfo));
  linkinfo->no=0;
  linkinfo->parent=NULL;
  linkinfo->next=NULL;
  linkinfo->prev=NULL;

  ((rsPtr)currentNode->_private)->linkinfo=linkinfo;

  return (0);
}
  

int createprivate(xmlNodePtr currentNode){

  rsPtr private=NULL;

  private=(rsPtr)malloc(sizeof(rels));
  private->linksPtr=NULL;
  private->linksOrigin=NULL;
  private->associations=NULL;
  private->local=NULL;
  private->copies=NULL;
  private->linkinfo=NULL;
  private->mrelations=NULL;

  currentNode->_private=private;

  return(0);
}

int axunlink(xmlNodePtr currentNode){

  int i=0;
  rsPtr private=NULL;

  if(!currentNode)return(-1);
  if(!currentNode->_private)return(0);
  if(!((rsPtr)currentNode->_private)->linkinfo)return(0);

  for(i=0;i<((rsPtr)currentNode->_private)->linkinfo->no;i++){
    private=(rsPtr)currentNode->last->_private;

    /* link */

    if(private->linkinfo->prev)private->linkinfo->prev->next=currentNode->last;
    else private->linkinfo->parent->children=currentNode->last;

    if(private->linkinfo->next)private->linkinfo->next->prev=currentNode->last;

    /* unlink */

    currentNode->last->parent=private->linkinfo->parent;
    currentNode->last->next=private->linkinfo->next;

    if(currentNode->last->prev){
      currentNode->last=currentNode->last->prev;
      currentNode->last->next->prev=private->linkinfo->prev;
      currentNode->last->next=NULL;
    }
    else{
      currentNode->last->prev=private->linkinfo->prev;
      currentNode->last=NULL;
      currentNode->children=NULL;
    }

    free(private->linkinfo);
    private->linkinfo=NULL;
    
  }

  free(((rsPtr)currentNode->_private)->linkinfo);
  ((rsPtr)currentNode->_private)->linkinfo=NULL;

  return(0);
}

int linkup(xmlNodePtr currentNode){
  
  int i=1,retval=0,start=-1,end=-1;
  xmlNodePtr node=NULL;
  struct xpathObjWrap *linkNodes=NULL;
  
  if(!currentNode)return(0);
  if(!currentNode->_private)return(0);

  retval=getLinks(currentNode,&linkNodes);

  if(!linkNodes)return(0);

  locateNode(linkNodes, i, &node, &start, &end);

  for( i = 2; node; i++ ){

    if((start==-1)&&(end==-1)){

      /* link in node */

      if(!((rsPtr)currentNode->_private)->linkinfo)createlinkinfo(currentNode);

      ((rsPtr)currentNode->_private)->linkinfo->no++;
      
      /* store original node properties */

      if(!node->_private)createprivate(node);
      if(!((rsPtr)node->_private)->linkinfo)createlinkinfo(node);

      ((rsPtr)node->_private)->linkinfo->parent=node->parent;
      ((rsPtr)node->_private)->linkinfo->next=node->next;
      ((rsPtr)node->_private)->linkinfo->prev=node->prev;

      /* unlink */

      if(node->prev)node->prev->next=node->next;
      else(node->parent->children=node->next);
      if(node->next)node->next->prev=node->prev;

      /* link */

      node->parent=currentNode;
      node->prev=currentNode->last;
      node->next=NULL;
      if(node->prev)node->prev->next=node;
      else currentNode->children=node;
      currentNode->last=node;

    }
    locateNode(linkNodes, i, &node, &start, &end);
  }

  freeWrap (&linkNodes);
  return(0);
}
      

int getLinks(xmlNodePtr currentNode, struct xpathObjWrap** linkNodes){

  lsPtr linksLocal=NULL;
  struct xml_data* xml_data_local=NULL;
  xmlNodePtr localNode=NULL;
  xmlNodePtr node=NULL;
  struct xpathObjWrap* keys=NULL;
  int i=0;

  if((!xmlStrcmp(currentNode->name,(xmlChar*)"complex"))&&(!xmlStrcmp(currentNode->ns->href,(xmlChar*)"http://www.grids.ac.uk/eccp/ns#"))){

    if(currentNode->_private){
      if((linksLocal=((rsPtr)currentNode->_private)->linksOrigin)){
	while(linksLocal){
	  linkWrap(linkNodes,linksLocal->nodes);
	  linksLocal=linksLocal->next;
	}
      }
      if(((rsPtr)currentNode->_private)->local){
	linkWrap(linkNodes,(((rsPtr)currentNode->_private)->local));
      }
      if(((rsPtr)currentNode->_private)->associations){
	linkWrap(linkNodes,(((rsPtr)currentNode->_private)->associations));
      }
    }
  }
  
  else{
    
    if(currentNode->_private){
      
      xml_data_local=gdata->prev;
      
      while(xml_data_local){
	locateNode(*(xml_data_local->set),*(xml_data_local->set_pos),&localNode,(NULL),(NULL));
	if(localNode->_private){
	  if(((rsPtr)localNode->_private)->local){
	    getKeys(&keys,((rsPtr)localNode->_private)->local);
	    if(keys){

	      i = 1;
	      locateNode(keys, i, &node, (NULL), (NULL));

	      for( i = 2; node; i++){

		linksLocal=((rsPtr)currentNode->_private)->linksOrigin;
		while(linksLocal){
		  if((linksLocal->key==node)||(linksLocal->scope==2)){
		    if(((linksLocal->start==-1)&&(linksLocal->end==-1))||
		       ((linksLocal->start<=gdata->sel_index)&&(linksLocal->end>=gdata->sel_index2))){
		      linkWrap(linkNodes,linksLocal->nodes);
		    }
		  }
		  linksLocal=linksLocal->next;
		}
		locateNode(keys, i, &node, (NULL), (NULL));
	      }
	    }
	    break;
	  }
	}
	xml_data_local=xml_data_local->prev;
      }
      
      if(!xml_data_local){
	linksLocal=((rsPtr)currentNode->_private)->linksOrigin;
	while(linksLocal){
	  if(linksLocal->scope!=1){
	    linkWrap(linkNodes,linksLocal->nodes);
	  }
	  linksLocal=linksLocal->next;
	}
      }
      linkWrap(linkNodes,(((rsPtr)currentNode->_private)->associations)); 
    }
  }
  return(0);
}


int getMLinks(){
  return(0);
}


xmlXPathObjectPtr find(xmlChar* resource,xmlChar* xptrExpr){

  xmlNodePtr targetNode=NULL,linkNode=NULL;
  struct xpathObjWrap* targetNodes=NULL;
  xmlXPathObjectPtr matchingObject=NULL;
  xmlNodeSetPtr matchingNodeSet=NULL;
  int i=0,j=0,k=0;
  xmlNodePtr context=NULL;
  xmlNodePtr currentNode=NULL;
  struct xpathObjWrap* linkNodes=NULL;
  int retval=0;
  struct ctxStr* ctxStrLocal=NULL;

  if(!xptrExpr)return(NULL);
  if(!resource)return(NULL);
  if(!gdata)return(NULL);

  context=NULL;
  ctxStrLocal=ctxStrOrigin;
  while(ctxStrLocal){
    if(ctxStrLocal->type==1){
      if(xptr_query(ctxStrLocal,ctxStrLocal->original_doc_root,xptrExpr))return(NULL);
      linkWrap(&targetNodes,ctxStrLocal->xpathObj);
      freeWrap(&ctxStrLocal->xpathObj);
      ctxStrLocal->xpathObj=NULL;
    }
    ctxStrLocal=ctxStrLocal->next;
  }
  
  if(!targetNodes)goto cleanup;
  
  if(!linksResolved)resolveLinks();

  currentNode=gdata->sel;
  
  /* getLinks and getMLinks determine matching nodes */
  /* linkNodes contains only those matching nodes */
  
  getLinks(currentNode,&linkNodes);

  getMLinks(currentNode,&linkNodes);

  /* process the links structure and find matching nodes */
    
    if(linkNodes){
      if(!matchingNodeSet){
	matchingNodeSet=(xmlNodeSetPtr)malloc(sizeof(xmlNodeSet));
	matchingNodeSet->nodeTab=NULL;
      }

      i = 1;
      locateNode(targetNodes, i, &targetNode, (NULL), (NULL));

      for(i = 2; targetNode; i++){

	j = 1;
	locateNode(linkNodes, j, &linkNode, (NULL), (NULL));

	for(j = 2; linkNode; j++){

	  if(targetNode==linkNode){
	    matchingNodeSet->nodeTab=realloc(matchingNodeSet->nodeTab,(k+2)*sizeof(xmlNodePtr));
	    matchingNodeSet->nodeTab[k]=targetNode;
	    k++;
	  }
	  locateNode(linkNodes, j, &linkNode, (NULL), (NULL));
	}
	locateNode(targetNodes, i, &targetNode, (NULL), (NULL));
      }
    }
  
    if(k!=0){

      matchingObject=(xmlXPathObjectPtr)malloc(sizeof(xmlXPathObject));     
      matchingObject->nodesetval=NULL;
      matchingObject->type=1;
      matchingObject->boolval=0;
      matchingObject->nodesetval=matchingNodeSet;
      matchingObject->nodesetval->nodeNr=k;
      matchingObject->nodesetval->nodeMax=k;
    }
    
  cleanup:

    if (targetNodes) freeWrap (&targetNodes);
    if (linkNodes) freeWrap (&linkNodes);
    if (retval == -1) k = retval;
    
    return matchingObject;
    
}


int noTracerProperty(){

  struct xpathObjWrap* objWrap = NULL;
  xmlChar **card = NULL, **cardOrig = NULL, **property = NULL, **propertyOrig = NULL, *item = NULL;
  int cardValue = 0, retval = 0, no = 0;

  /* find a property with a non zero cardinality */

  startResource(&objWrap,gdata->class);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2000/01/rdf-schema#subClassOf",0);
  RDFLiteral(&cardOrig,(xmlChar*)"http://www.w3.org/2002/07/owl#cardinality",objWrap);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2002/07/owl#cardinality",0);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2002/07/owl#cardinality",1);
  RDFTrace(&objWrap,(xmlChar*)"http://www.w3.org/2002/07/owl#onProperty",0);
  RDFResource(&propertyOrig,objWrap);

  if ((propertyOrig) && (cardOrig)){
    property = propertyOrig;
    card = cardOrig;
    while (card){
      if (card && (**card != '0')) break;
      card = card + 1;
      property = property + 1;
    }
    if (card) sscanf((char*)*card, "%d", &cardValue);
    else ERROR;
  }

  else ERROR;

  if (property) item = *property + xmlStrlen(base_uri_string);
  else ERROR;

  if (cardValue > 0){
    while((no = xml_get(item)) > 0){
      retval += no;
      axDeselect();
      gdata->data_sets_position++;
    }
    if(retval > 0) retval /= cardValue;
  }
  
 cleanup:
  
  if (propertyOrig) freeStringArray(&propertyOrig);
  if (cardOrig) freeStringArray(&cardOrig);
  freeWrap (&objWrap);

  return(retval);
}


int ctx_xml_select(struct ctxStr* ctxStrLocal,xmlChar* term,int type){

  int retval=0;
  xmlChar* resource=NULL;

  if(!linksResolved)resolveLinks();

  if(ldata){
    if(ctxStrLocal->xml_data_ptr->next_resource){
      free(ctxStrLocal->xml_data_ptr->next_resource);
    }
    ctxStrLocal->xml_data_ptr->next_resource=(xmlChar*)malloc((xmlStrlen(base_uri_string)+xmlStrlen(term)+1)*sizeof(xmlChar));
    resource=ctxStrLocal->xml_data_ptr->next_resource;
  }
  
  else{
    if(ctxStrLocal->next_resource){
      free(ctxStrLocal->next_resource);
    }
    ctxStrLocal->next_resource=(xmlChar*)malloc((xmlStrlen(base_uri_string)+xmlStrlen(term)+1)*sizeof(xmlChar));
    resource=ctxStrLocal->next_resource;
  }
  
  sprintf((char*)resource,"%s%s",(char*)base_uri_string,(char*)term);

  /* what are we selecting? */

  if (!type) type = findType();

  switch (type){
  case(0):
    /* bail out if no information found */
    retval=-1;
    break;
  case(1):
    /* handle Class selection */
    retval=xml_select_ds(ctxStrLocal,term);
    break;
  case(2):
    /* handle DatatypeProperty selection */
    retval=xml_get(term);
    break;
  case(3):
    /* handle ObjectProperty selection */
    retval=xml_select_ds(ctxStrLocal,term);
    break;
  case(4):
    /* handle FunctionalProperty selection */
    retval=xml_get(term);
    break;
  }

  return retval;
}


int xml_select_ds(struct ctxStr* ctxStrLocal,xmlChar* term){

  xmlXPathObjectPtr result = NULL;
  xmlChar *xptrExpr = NULL, **xlocator = NULL, **resource = NULL;
  xmlChar **dataSetType = NULL, **delimiter = NULL;
  int retval = 0, len = 0, dilink = 0, dstype = 0, found = 0;
  xmlNodePtr context = NULL;
  struct xpathObjWrap *evaluation = NULL;
  struct classStr* classStrPtr = NULL;

  if (ldata) resource = &(ctxStrLocal->xml_data_ptr->next_resource);
  else resource = &(ctxStrLocal->next_resource);

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));

  getClassMap(&classStrPtr, *resource);

  if ( classStrPtr ){
    
    xlocator = classStrPtr->xlocator;
    delimiter = classStrPtr->delimiter;
    dstype = classStrPtr->dataSetType;
    
  }

  if(gdata){
    if(gdata->sel){
      if(gdata->sel->type!=XML_ELEMENT_NODE){
	dilink=1;
      }
    }
  }
    
  if(!dilink){
    context=ctxStrLocal->original_doc_root;
      
    if(ldata){
      if(ldata->sel){
	context=ldata->sel;	  
      }
    }
      
    len=xmlStrlen(*resource)+41;
    if(len>99){
      xptrExpr=(xmlChar*)realloc(xptrExpr,(len+1)*sizeof(xmlChar));
    }
    
    sprintf((char*)xptrExpr,"xpointer(//eccp:complex[@xlink:role='%s'])",(char*)*resource);
    xptr_query(ctxStrLocal,context,xptrExpr);
    linkWrap(&evaluation,ctxStrLocal->xpathObj);
    freeWrap(&ctxStrLocal->xpathObj);
    ctxStrLocal->xpathObj=NULL;
    
  }

  if(xlocator){
    
    len=xmlStrlen(*xlocator);
    if(len>99){
      xptrExpr=(xmlChar*)realloc(xptrExpr,(len+1)*sizeof(xmlChar));
    }
    
    sprintf((char*)xptrExpr,"%s",(char*)*xlocator);
    if(!dilink){
      xptr_query(ctxStrLocal,context,xptrExpr);
      linkWrap(&evaluation,ctxStrLocal->xpathObj);
      freeWrap(&ctxStrLocal->xpathObj);
      ctxStrLocal->xpathObj=NULL;
    }
  }
  
  if(dilink){
    if(ldata){
      if ( (result = find(*resource, xptrExpr)) ) wrapObject(&evaluation, result);
    }
  }
  
  if (evaluation){
    
    create_xml_data(ctxStrLocal);
    
    if (ldata->class) free(ldata->class);
    ldata->class = *resource;
    if (xlocator)ldata->xlocator = *xlocator;
    ldata->set = &ldata->data_sets;
    ldata->set_pos = &ldata->data_sets_position;
    ldata->data_sets = evaluation;
    ldata->delimiter = delimiter;
    locateNode(ldata->data_sets, ldata->data_sets_position, &ldata->sel, (NULL), (NULL));
    
    retval = countNodes(ldata->data_sets);

    /* find number of logical data sets */

    if (dstype == 1){
      gdata->data_sets_type = dstype;
      retval = noTracerProperty();
      if (retval == -1) {
	axDeselect();
	goto cleanup;
      }
      ldata->data_sets_position = 1;
    }
    
    ldata->data_sets_count = retval;

  }
    
  cleanup:
    
  if (resource) *resource = NULL;
  if (xptrExpr) free(xptrExpr);
  
  return retval;
  
}

int xml_select(xmlChar* term,int type){

  int retval=0,sum=0;
  struct ctxStr* ctxStrLocal=NULL;

  if(!ctxStrOrigin)return(-1);
  if(!ctxStrPtr)return(-1);
  
  if(!gdata){

    ctxStrPtr = ctxStrOrigin;
    while(ctxStrPtr){
      if(ctxStrPtr->type == 1){
	retval = ctx_xml_select(ctxStrPtr, term, type);
	if (retval == -1) return(-1);
	else{
	  sum += retval;
	  if (retval) {
	    if (!ctxStrLocal) ctxStrLocal = ctxStrPtr;
	    else *(ctxStrPtr->xml_data_ptr->set_pos) = 0;
	  }
	}
      }
      ctxStrPtr = ctxStrPtr->next;
    }
    if (sum == 0) ctxStrPtr = ctxStrOrigin;
    else ctxStrPtr = ctxStrLocal;
  }
  
  else {sum = ctx_xml_select(ctxStrPtr, term, type);}
  
  return sum;

}


xmlNodePtr resolveNode(struct ctxStr* ctxStrLocal,xmlNodePtr copyNode){

  xmlNodePtr nodeCopy=NULL;
  xmlChar *resVars=NULL;
  xmlChar *resExpr=NULL;
  xmlAttrPtr id=NULL;
  xmlChar *xptrExpr=NULL;
  xmlAttrPtr propLocal=NULL;

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));      

  /* make a copy of the node */
  
  nodeCopy=xmlCopyNode(copyNode,2);
  nodeCopy->doc=copyNode->doc;
  nodeCopy->parent=copyNode->parent;
  nodeCopy->next=copyNode->next;
  
  /* link node into main tree */
  
  nodeCopy->next->prev=nodeCopy;
  
  /* get id attributes */
  
  propLocal=nodeCopy->properties;
  while(propLocal){
    if(propLocal->name){
      if(!xmlStrcmp(propLocal->name,(xmlChar*)"id")){
	break;
      }
    }
    propLocal=propLocal->next;
  }
  
  id=propLocal;
  
  /* resolve id attribute value */
  
  if(id){
    resolveVars(ctxStrLocal,id->children->content,&resVars);
    resolveExpr(resVars,&resExpr);
    free(id->children->content);
    id->children->content=resExpr;
    free(resVars);
    resVars=NULL;
  }

  if(xptrExpr)free(xptrExpr);
  return nodeCopy;
}


int expandDOM(struct ctxStr* ctxStrLocal){

  /* all DOMs should be expanded first so that evaluations can be made across documents */
  /* evaluations are left to resolveLink() */
  /* link and complex elements should be expanded */

  xmlNodePtr copyNode,copyNodePrev=NULL;
  int linkinstances=0,retval=0,i=1,j=0;
  xmlChar* inst=NULL;
  xmlNodePtr nodeCopy=NULL,context=NULL;
  xmlChar *xptrExpr=NULL;
  struct xpathObjWrap *copyNodes=NULL;

  if(!ctxStrLocal) ERROR;

  xptrExpr=(xmlChar*)malloc(100*sizeof(xmlChar));

  /* find link and complex nodes */

  context = ctxStrLocal->original_doc_root;
  sprintf((char*)xptrExpr, "xpointer(//eccp:complex | //eccp:link | //molpro:association)");
  if (xptr_query(ctxStrLocal,context,xptrExpr)) ERROR;
  linkWrap(&copyNodes,ctxStrLocal->xpathObj);
  freeWrap(&ctxStrLocal->xpathObj);
  ctxStrLocal->xpathObj = NULL;
  linkWrap(&ctxStrLocal->copyNodes,copyNodes);
  
  if(copyNodes){

    i = 1;
    locateNode(copyNodes, i, &copyNode, (NULL), (NULL));

    for( i = 2; copyNode; i++){

      copyNodePrev=copyNode->prev;

      /* get value of the instances attribute */

      inst=getProperty(copyNode,(xmlChar*)"instances",(NULL));
      
      /* if instances attribute not found, default to 1 */

      if(!inst)linkinstances=1;
      else {
	sscanf((char*)inst,"%d",&linkinstances);
	inst=NULL;
      }

      /* make copies of the node and store addresses */

      for(j=1;j<=linkinstances;j++){

	if(!ctxStrLocal->var)ctxStrLocal->var=(xmlChar*)malloc(100*sizeof(char));
	sprintf((char*)ctxStrLocal->var,"%d",j);
	  
	nodeCopy=resolveNode(ctxStrLocal,copyNode);
	  
	/* complete tree linking */
	
	nodeCopy->prev=copyNodePrev;
	nodeCopy->prev->next=nodeCopy;
	if(nodeCopy->prev==NULL)nodeCopy->parent->children=nodeCopy;
	
	copyNodePrev=nodeCopy;

	nodeCopy->_private=NULL;

	if(!copyNode->_private){
	  copyNode->_private=(rsPtr)malloc(sizeof(rels));
	  ((rsPtr)copyNode->_private)->copies=NULL;
	}
	((rsPtr)copyNode->_private)->copies=realloc(((rsPtr)copyNode->_private)->copies,(j+1)*sizeof(xmlNodePtr));
	((rsPtr)copyNode->_private)->copies[j-1]=nodeCopy;
	((rsPtr)copyNode->_private)->copies[j]=NULL;
      }
      locateNode(copyNodes, i, &copyNode, (NULL), (NULL));
    }
  }
    
 cleanup:  
    
  freeWrap(&copyNodes);
  if(xptrExpr)free(xptrExpr);
  
  return retval;
}


int uri_get(xmlChar* uri){

  xmlChar* buffer=NULL;
  xmlChar* filename=NULL;
  xmlDocPtr doc=NULL;
  struct ctxStr* ctxStrLocal=NULL;
  int retval=0;
  void* httpContext=NULL;

  if(!uri)return(-1);

  xpathCtx=NULL;

  if(!pctxt){
    pctxt = xmlNewParserCtxt();
    if (!pctxt) ERROR;
    /*    xmlCtxtUseOptions(pctxt,0);*/
  }
  
  doc=xmlCtxtReadFile(pctxt,(char*)uri,NULL,0);

  /* Create_xpath evaluation context */
  
  if (!doc) ERROR;

  xpathCtx = xmlXPtrNewContext(doc,NULL,NULL);
  if(!xpathCtx) ERROR;

  ctxStrLocal=ctxStrOrigin;

  if(!ctxStrLocal){
    ctxStrOrigin=(struct ctxStr*)malloc(sizeof(struct ctxStr));
    ctxStrLocal=ctxStrOrigin;
    ctxStrPtr=ctxStrOrigin;
    ctxStrLocal->prev=NULL;}
    
  else{
    while(ctxStrLocal->next){
      ctxStrLocal=ctxStrLocal->next;
    }
    ctxStrLocal->next=(struct ctxStr*)malloc(sizeof(struct ctxStr));
    ctxStrLocal->next->prev=ctxStrLocal;
    ctxStrLocal=ctxStrLocal->next;
  }
    
  ctxStrLocal->next=NULL;
  ctxStrLocal->xpathObj=NULL;
  ctxStrLocal->xpathCtx=xpathCtx;
  ctxStrLocal->parent=NULL;
  ldata=NULL;
  ctxStrLocal->xml_data_origin=NULL;
  ctxStrLocal->original_doc_root=xmlDocGetRootElement(xpathCtx->doc);
  ctxStrLocal->next_resource=NULL;
  ctxStrLocal->var=NULL;
  ctxStrLocal->var2=NULL;
  ctxStrLocal->copyNodes=NULL;

 cleanup:

  if(filename)free(filename);
  if(buffer)free(buffer);
  if(retval==-1){
    if(doc){
      xmlFreeDoc(doc);
      doc=NULL;
    }
    if(xpathCtx){
      xmlXPathFreeContext(xpathCtx);
      xpathCtx=NULL;
    }
  }
  if(httpContext)xmlNanoHTTPClose(httpContext);

  return(retval);
}

int loadDoc(){

  xmlChar* xptrExpr=NULL;
  int i=0;
  int docs_found=0;
  xmlChar** uris=NULL;
  struct ctxStr* ctxStrLocal=NULL;
  struct ctxStr* parent=NULL;
  xmlNodePtr context=NULL;
  int retval=0;
  struct xpathObjWrap* results=NULL;

  xptrExpr=(xmlChar*)malloc(100*sizeof(char));

  ctxStrLocal=ctxStrOrigin;
  while(ctxStrLocal->next){
    ctxStrLocal=ctxStrLocal->next;
  }

  parent=ctxStrLocal;

  sprintf((char*)xptrExpr,"xpointer(//eccp:semantics/@location)");
  context=ctxStrLocal->original_doc_root;

  retval=type_query(4,context,xptrExpr);
  retval=collectResults(&results);
  retval=getDataElements(results,&uris,&docs_found,NULL);

  if(uris){
    supdocs=1;
    for(i=0;*(uris+i);i++){
      if(uri_get(*(uris+i))) continue;
      register_namespaces(xpathCtx);
      ctxStrLocal=ctxStrLocal->next;
      ctxStrLocal->type=2;
      ctxStrLocal->parent=parent;
      expandDOM(ctxStrLocal);
    }
  }

  if (results) freeWrap(&results);
  if (xptrExpr) free(xptrExpr);
  if (uris) free(uris);

  return (retval);
}


int axProxy(char* proxy){
  xmlChar* proxyuri=NULL;
  proxyuri=(xmlChar*)xmlStrdup((xmlChar*)proxy);
  return(0);
}


int xml_get_uri(xmlChar* uri){

  struct ctxStr* ctxStrLocal=NULL;

#ifdef DEBUG
  printf( "xml_get_uri: %s\n", uri );
#endif

  if(uri_get(uri))return(-1);
  register_namespaces(xpathCtx);
  
  ctxStrLocal=ctxStrOrigin;
  while(ctxStrLocal->next){
    ctxStrLocal=ctxStrLocal->next;
  }

  ctxStrLocal->type=1;
  expandDOM(ctxStrLocal);

  loadDoc();

  return (0);
}


int axGetUri(char* uri){

  struct ctxStr* ctxStrLocal=NULL;

  if(uri_get((xmlChar*)uri))return(-1);
  register_namespaces(xpathCtx);

  ctxStrLocal=ctxStrOrigin;
  while(ctxStrLocal->next){
    ctxStrLocal=ctxStrLocal->next;
  }

  ctxStrLocal->type=3;
  expandDOM(ctxStrLocal);
  supdocs=1;

  return(0);

}


void register_namespaces(xmlXPathContextPtr xpathCtx){

  /* register a standard set of namespaces
     and internal functions */
  
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"xlink",(xmlChar*)"http://www.w3.org/1999/xlink");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"xsd",(xmlChar*)"http://www.w3.org/2001/XMLSchema#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"dc",(xmlChar*)"http://purl.org/dc/elements/1.1/");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"rdfs",(xmlChar*)"http://www.w3.org/2000/01/rdf-schema#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"daml",(xmlChar*)"http://www.daml.org/2001/03/daml+oil#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"rdf",(xmlChar*)"http://www.w3.org/1999/02/22-rdf-syntax-ns#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"owl",(xmlChar*)"http://www.w3.org/2002/07/owl#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"eccp",(xmlChar*)"http://www.grids.ac.uk/eccp/ns#");
  xmlXPathRegisterNs(xpathCtx,(xmlChar*)"molpro",(xmlChar*)"http://www.molpro.net/schema/molpro2005");

  registerInternalFunctions(xpathCtx);

}

int axCount(){

  struct xml_data *dptr = NULL;
  int count = 0;

  return ( ( ctxStrPtr && ( dptr = ctxStrPtr->xml_data_ptr ) ) ? dptr->data_sets_count : 0 );
}

/* return the number of nodes in the logical path */

int axDepth(){

  int depth = 0;
  struct xml_data *dptr = ctxStrPtr->xml_data_origin;

  while( dptr ){
    depth++;
    if( dptr->set == &(dptr->item_set) ) depth++;
    dptr = dptr->next;
  }

  return depth;
}

/* reset the query context to 'document' */

int axReset()
{

  int depth = 0;
  
  depth = axDepth();
  while( depth-- ) axDeselect();
  
  return 0;
}

/* return the position of the selected data set */

int axPosition(){

  struct xml_data *dptr = NULL;

  return ( ( ctxStrPtr && ( dptr = ctxStrPtr->xml_data_ptr ) ) ? dptr->data_sets_position : 0 );
  
}

