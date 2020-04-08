
/* AgentX error functions */

#include<string.h>
#include<stdio.h>
#include<axError.h>
#include<axArray.h>

static axArray_t *errorBuff = NULL;
static FILE *errorFile = NULL;

int axCheckBufLen( char *str, unsigned int size ){

  if ( ! str ) return -1;

  if ( strlen( str ) > size ) {
    axRegError( "AgentX error: Buffer overflow\n" );
    return -1;
  }

  return 0;

}

/* Add to an array of error messages */

int axRegError( char *error ){

  return axDupToStringArray( &errorBuff, error );

}

/* Free error message memory */

int axClearErrors(){

  axFreeStringArray( errorBuff );
  errorBuff = NULL;

  return 0;
  
}


/* write errors to error stream */

int axPrintErrors(){

  int i, imax;
  char **stringArray = NULL;

  if ( ! errorBuff ) return 0;
  if ( ! errorFile ) errorFile = ERROR_SP;
  
  stringArray = (char**)axArrayData( errorBuff );
  imax = axArraySize( errorBuff );
  
  for ( i = 0; i < imax; i++ ) fprintf( errorFile, "%s\n", stringArray[ i ] );

  return 0;
}

/* Set error stream */

int axSetErrorStream( FILE *es ){

  if ( ! es ) return -1;
  errorFile = es;
  
  return 0;
}

/* get all errors */

int axGetErrors( char ***errors ){

  if( ! errors ) return -1;
  
  *errors = (char**)axArrayData( errorBuff );

  return ( axArraySize( errorBuff ) );
}

/* get top level error */

int axGetError( char **error ){

  char **arrayData = NULL;

  if ( ! error ) return -1;

  arrayData = axArrayData( errorBuff );
  
  if( arrayData ) *error = ( (char**)axArrayData( errorBuff ) )[ 0 ];
  
  return 0;
}
