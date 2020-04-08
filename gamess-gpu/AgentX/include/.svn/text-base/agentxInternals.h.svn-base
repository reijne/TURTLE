/* Abandon hope, ye who enters here... */
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <libxml/parser.h>
#include <libxml/xpointer.h>
#include <libxml/xpath.h>
#include <libxml/nanohttp.h>
#include <libxml/xpathInternals.h>

/* specify structures */

struct xpathObjWrap{
  int refs;
  xmlXPathObjectPtr xpathObj;
  struct xpathObjWrap **wrapTab;
  int tabNr;
};

struct classStr{
  xmlChar **xlocator;
  xmlChar **delimiter;
  int dataSetType;
};

struct propertyStr{
  xmlChar **xlocator;
  xmlChar **delimiter;
  xmlChar **sd;
  xmlChar **startIndex;
  xmlChar **endIndex;
  int startChar;
  int endChar;
};

struct evalCache{
  struct evalCache* next;
  struct evalCache* prev;
  xmlChar* expr1;
  xmlChar* expr2;
  void* load;
};

struct cacheBtree{
  int hash;
  struct evalCache* evalCacheOrigin;
  struct evalCache* evalCachePtr;
  struct cacheBtree* bigger;
  struct cacheBtree* smaller;
};

struct mLinksStr{
  int **to;
 struct mLinksStr* next;
  struct mLinksStr* prev;
  xmlNodePtr key;
  int scope;
  int start;
  int end;
};

struct mAssocStr{
  int *from;
  int **assoc;
  struct mAssocStr* next;
  struct mAssocStr* prev;
  int **local;
};

struct mRelationStr{
  struct mAssocStr* associations;
  struct mLinksStr* linksPtr;
  struct mLinksStr* linksOrigin;
  struct mRelationStr* next;
  struct mRelationStr* prev;
};
  
struct linksStr{
  struct linksStr* next;
  struct linksStr* prev;
  int start;
  int end;
  xmlNodePtr key;
  int scope;
  struct xpathObjWrap* nodes;
};

struct linkupinfo{
  int no;
  xmlNodePtr parent;
  xmlNodePtr next;
  xmlNodePtr prev;
};

struct relationStr{
  xmlNodePtr* copies;
  struct xpathObjWrap* local;
  struct xpathObjWrap* associations;
  struct linksStr* linksPtr;
  struct linksStr* linksOrigin;
  struct mRelationStr* mrelations;
  struct linkupinfo* linkinfo;
};

struct ctxStr{
  char* URI;
  int map_index;
  int persistent;
  struct mRelationStr* mrelations;
  xmlNodePtr original_doc_root;
  struct ctxStr* next;
  struct ctxStr* prev;
  struct ctxStr* parent;
  xmlXPathContextPtr xpathCtx;
  int type;
  struct xpathObjWrap *xpathObj;
  struct xml_data* xml_data_ptr;
  struct xml_data* xml_data_origin;
  xmlChar* next_resource;
  xmlChar* var;
  xmlChar* var2;
  struct xpathObjWrap *copyNodes;
};

/* some of this should be transfered to a parser context
 - work in progress */

struct xml_data{

  /* ontology information 
     URIs */

  xmlChar* class;
  xmlChar* property;

  /* variables */

  xmlChar* next_resource;
  xmlChar* xlocator;

  /* string pointers */

  xmlChar** cdata_list;
  xmlChar* charPtr;
  xmlChar* charData;
  xmlChar* bufferPtr;
  xmlChar* bufferData;
  xmlChar** freeable;

  /* data item information */

  struct xpathObjWrap *item_set;
  int item_set_position;
  int item_count;
  int item_pos;
  int card;
  int skip_delimiter;
  xmlChar** delimiter;
  xmlChar** startIndex;
  xmlChar** endIndex;
  int startChar;
  int endChar;
  int index;

  /* data set information */

  struct xpathObjWrap *data_sets;
  int data_sets_position;
  int data_sets_count;
  int data_sets_type;

  /* set points to either data_sets or item_set
     set_pos points to either data_sets_position or 
     item_set_position */

  struct xpathObjWrap **set;
  int *set_pos;

  /* sel points to the currently selected node
     the indices provide the location indices */

  xmlNodePtr sel;
  int sel_index;
  int sel_index2;

  struct xml_data *prev;
  struct xml_data *next;
};

#define ldata ctxStrLocal->xml_data_ptr
#define gdata ctxStrPtr->xml_data_ptr
#define ERROR {retval = -1; goto cleanup;}
#define END {retval = 0; goto cleanup;}

#define AX_INIT 1
#define AX_CLASS 2
#define AX_PROP 4
#define AX_SELECTION 8

#define AX_ON 1
#define AX_OFF 0

#define MAX_NS_LEN 100

/* Macros relating to location of INI file */
#define MAX_CONFIG_FNAME_LENGTH     512
#define CONFIG_FILENAME             "agentx.ini"
#define CONFIG_DIRECTORY            ".agentx"

typedef struct relationStr rels;
typedef struct relationStr* rsPtr;
typedef struct linksStr links;
typedef struct linksStr* lsPtr;

int axBuffer(int);
char* axBufferValue();
int axSelectType(char*,int);
int axNamespace(char*);
int axProxy(char*);
int axLen();  
int searchCache(struct cacheBtree* tree, struct cacheBtree** finalBranch, xmlChar* expr1, xmlChar* expr2, void** load, int* found);
int addToCache(struct cacheBtree** tree, xmlChar* expr1, xmlChar* expr2, void* load);
int ax_strstr(xmlChar *start, xmlChar *str, xmlChar **delimiter, int *sc, xmlChar** begin, xmlChar** end, int a, int b);
int xml_get_uri(xmlChar* uri);   /* select an XML data document */
int xml_select(xmlChar* term,int type);   /* select a dataset relating to a concept */
int xml_get(xmlChar* property);   /* return a data item */
int xml_get_next(int dir);   /* return next data item */
int xml_select_next(int dir);   /* select next dataset */
int xml_deselect();   /* drop the previous dataset selection */
int xml_refine(const xmlChar* expression);   /* refine a previous selection */

void register_namespaces(xmlXPathContextPtr xpathCtx);   /* register a standard set of namespaces with xpath contexts */
int uri_get(xmlChar* uri);   /* retrieve a document */
int cdata(int type,xmlNodePtr context,xmlChar* xptrExpr,xmlChar*** cdata_list,int* cdata_list_index,xmlChar*** freeable);   /* return the character data associated with the evaluation of the XPointer expression xptrExpr */
int loadDoc();   /* load RDF documents requested in XML documents */
xmlChar* locationCastToString(xmlXPathObjectPtr location);   /* return character data associated with the XPointer location */
xmlXPathObjectPtr find(xmlChar* resource,xmlChar* xptrExpr);   /* find a dataset related, to a current dataset selection, with XLink */
void message_handler(xmlChar* function_name,xmlChar* message,int type);   /* error handling routine */
int resolveLinks();
float evaluateExpression(char* expression);
int resolveExpr(xmlChar* expr,xmlChar** resExpr);
int addXPathObjects(xmlXPathObjectPtr* sum,xmlXPathObjectPtr object);

int destroy_xml_data(struct ctxStr* ctxStrLocal);
int create_xml_data(struct ctxStr* ctxStrLocal);
int xptr_query(struct ctxStr* ctxStrLocal,xmlNodePtr context,xmlChar* xptrExpr);
int resolveXPtrExpr(struct ctxStr* ctxStrLocal,xmlChar** xptrExpr,xmlChar** resExpr);
int ctx_xml_select(struct ctxStr* ctxStrLocal,xmlChar* term,int type);

int axLocateNode(struct xpathObjWrap* objWrap, int *no, xmlNodePtr* node, int* start, int* end);
int locateNode(struct xpathObjWrap* objWrap, int no, xmlNodePtr* node, int* start, int* end);

int initItems();
xmlChar* getProperty(xmlNodePtr node,xmlChar* property,xmlChar* ns);
int getNodes(int inst, xmlChar* exp, struct ctxStr* ctxStrLocal, struct xpathObjWrap** nodes);
int linkup(xmlNodePtr currentNode);
int getevalcontext(xmlNodePtr context,xmlChar** evalContext);
int dutilunlink(xmlNodePtr currentNode);
int createprivate(xmlNodePtr currentNode);
int getLinks(xmlNodePtr currentNode, struct xpathObjWrap **linkNodes);
int xml_select_ds(struct ctxStr* ctxStrLocal,xmlChar* term);
int axunlink(xmlNodePtr currentNode);
int checkPropertyRequest();
int axCheckState(int, int);
int count();
int depth();
int reset();
int position();
       
