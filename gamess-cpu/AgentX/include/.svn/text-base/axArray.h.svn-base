/* AgentX generic array handling functions */

#define BUF_INC 10

#ifndef axarraydef

typedef struct _axArray{
  unsigned int size;
  void *data;
} axArray_t;

#define axarraydef

#endif

int axFreeStringArray( axArray_t *axArray );
int axFreePtrArray( axArray_t *axArray );
void *axArrayData( axArray_t *axArray );
int axArraySize( axArray_t *axArray );
int axDupToStringArray( axArray_t **axArray, char *nstr );
int axChkPtrArray( axArray_t **axArray );
int axNewArray( axArray_t **axArray );
int axChkArray( axArray_t **axArray, int bytes );
int axLinkToStringArray( axArray_t **axArray, char *nstr );
int axDupToIntArray( axArray_t **axArray, int value );
void *axPtrArrayElement( axArray_t *axArray, int no );
int axSetPtrArrayElement( axArray_t **axArray, int no, void* nptr );
void *axGetPtrArrayElement( axArray_t *axArray, int no );
int axInitPtrArray( axArray_t **axArray, unsigned int size );
