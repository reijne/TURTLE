
/* AgentX generic array handling functions */

#include<axArray.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

/* array creation/ extension functions */

int axDupToIntArray( axArray_t **axArray, int value ){
  
  int size;
  int *data;
  
  if ( axChkArray( axArray, sizeof( int ) ) < 0 ) return -1;
  
  size = axArraySize( *axArray );
  data = (int*)axArrayData( *axArray );
  data[ size - 1 ] = value;

  return 0;
}

int axDupToStringArray( axArray_t **axArray, char *nstr ){
  
  int size;
  char **data;

  if ( axChkArray( axArray, sizeof( char* ) ) < 0 ) return -1;

  size = axArraySize( *axArray );
  data = (char**)axArrayData( *axArray );
  data[ size - 1 ] = strdup ( nstr );

  return 0;
}

int axLinkToPtrArray( axArray_t **axArray, void *nstr ){

  int size;
  void **data;
  
  if ( axChkArray( axArray, sizeof( void* ) ) < 0 ) return -1;

  size = axArraySize( *axArray );
  data = (void**)axArrayData( *axArray );
  data[ size - 1 ] = nstr;

  return 0;
}

/* array freeing functions */

int axFreeIntArray( axArray_t *axArray ){

  void *data;
  
  if( axArray && ( data = axArray->data ) ) free( data );
  
  return 0;

}

int axFreeStringArray( axArray_t *axArray ){

  int size;

  if ( ! axArray ) return -1;
  
  size = axArray->size;
  while( size-- ) free( ((char**)(axArray->data))[ size ] );
  free( axArray->data );
  
  free ( axArray );

  return 0;
}

int axFreePtrArray( axArray_t *axArray ){
  
  if( axArray ){
    free( axArray->data );
    free( axArray );
  }

  return 0;
}

/* create an initialised ptr array with size elements */

int axInitPtrArray( axArray_t **axArray, unsigned int size ){

  int allocSize;

  allocSize = (int)floor( (double)size/BUF_INC );

  axNewArray( axArray );
  (*axArray)->data = calloc( ( allocSize + 1 ) * BUF_INC, sizeof( void* ) );
  (*axArray)->size = size;

  return 0;
}


/* return the nth pointer of a pointer array */

void *axGetPtrArrayElement( axArray_t *axArray, int no ){

  if( ! axArray ) return NULL;
  if( no > axArraySize( axArray ) ) return NULL;

  return ( (void**)(axArray->data) )[ no - 1 ];

}
  
/* set the nth pointer of a pointer array ( allocate */
/* memory if required                                */

int axSetPtrArrayElement( axArray_t **axArray, int no, void* nptr ){

  int asize, inc;
  
  if( ! *axArray ) axNewArray( axArray );
  
  if( no > ( asize = axArraySize( *axArray ) ) ){
    inc = no - asize ;
    while( inc-- ) axLinkToPtrArray( axArray, NULL );
  }

  ( (void**)( axArrayData( *axArray ) ) )[ no - 1 ] = nptr;

  return 0;

}

/* return the array data */

void *axArrayData( axArray_t *axArray ){

  return ( axArray ? axArray->data : NULL );
}

/* return the number of array elements */

int axArraySize( axArray_t *axArray ){
  
  return ( axArray ? axArray->size : 0 );

}

/* allocate memory for a new array */

int axNewArray( axArray_t **axArray ){
  
  return ( ( *axArray = (axArray_t *)calloc( 1, sizeof( axArray_t ) ) ) ? 0 : -1 );

}

/* check memory allocation and array size */

int axChkArray( axArray_t **axArray, int bytes ){

  /* BUF_INC is set in axArray.h */

  unsigned int size;
  unsigned int *sizePtr = NULL;
  void **data = NULL;

  if ( ! *axArray ) axNewArray( axArray );
  data = &( (*axArray)->data );
  if ( ! *data ) *data = calloc( BUF_INC, bytes );

  sizePtr = &( (*axArray)->size );

  (*sizePtr) ++;
  size = *sizePtr;

  if ( ! ( size % BUF_INC ) ) *data = realloc( *data, ( size + BUF_INC ) * bytes );
  memset( (char*)(*data) + ( size * bytes ), 0, bytes ); 

  return 0;
}
