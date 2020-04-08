
/* AgentX logical path handling routines */

#include<axPath.h>
#include<agentx.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


int axSubstitute( char *path, char **str ){

  int retval = 0, i;
  char **strArray = NULL, *strptr = NULL;

  if( !str ) ERROR;
  if( !path ) ERROR;
  if( axPath( path, &strArray ) < 0 ) ERROR;

  if( strArray ){
    
    *str = *strArray;
    for( i = 1; strptr = strArray[ i ]; i++ ){
      free( strptr );
    }
    free( strArray );
  }
  
  else *str = NULL;  
  
 cleanup:
  
  return retval;
}


/* obtain the values of all properties in the context represented by path */

int axPath( char *path, char ***strArray ){
  
  int noValues = 0, retval = 0;
  struct _axArray *result = NULL;

  if( axPathNav( path, &result ) < 0 ) ERROR;
  *strArray = (char**)axArrayData( result );
  noValues = axArraySize( result );

 cleanup:
  axReset();
  if( result ) free( result );
  if( ! retval ) retval = noValues; 
  return retval;
} 


int axPathValue( char *prop, axArray_t **result ){

  char buf[30], *cp;
  int count, pos;

  /* support for 'count' and 'position' should be built into AgentX */
  
  if( strcmp( prop, "count" ) == 0 )
    {
      count = axCount();
      if( snprintf( buf, 30, "%d", count ) >= 30 )
	{
	  fprintf( stderr, "Buffer too small to hold data!\n" );
	  return -1;
	}
      axDupToStringArray( result, buf );
    }
  else if ( strcmp( prop, "position") == 0 )
    {
      pos = axPosition();
      if( snprintf( buf, 30, "%d", pos ) >= 30 )
	{
	  fprintf( stderr, "Buffer too small to hold data!\n" );
	  return -1;
	}
      axDupToStringArray( result, buf );
    }
  else
    {
      /* this assumes a single value - need to implement axValueAll() */
      if( ! ( cp = axValue() ) )
	{
	  fprintf( stderr, "Error selecting value!\n" );
	  return -1;
	}
      axDupToStringArray( result, cp );
    }

  return 0;
}


int axPathNav( char *path, struct _axArray **result )
{
  char *cp, buf[30];
  int pos, count, nEntity;  
  struct axOperation op; 

  cp = path;
  
#ifdef DEBUG
  printf( "axPath: property %s\n", path );
  printf( "axPath: initial context\n" );
  axCurrent();
#endif
  
  /* Walk logical path */
  while( *cp )
    {
      /* Parse next chunk of logical path and get corresponding axOperation */
      if( ( cp = convertToAxOperation( cp, &op ) ) == NULL ) return -1;
      
      if( *cp == ']' ) cp++;
      
      if( *cp == '.' ) cp++;

      /* Process axOperation */
      if( strcmp( op.entity, "count" ) && strcmp( op.entity, "position" ) )
        {
	  if( ( nEntity = doAxOperation( &op, cp, result ) ) < 0 ) return nEntity;
	}

      if( op.selection == -2 ) return nEntity;

      if( *cp == ',' )
	{
	  if( nEntity > 0 ) {
	    axPathValue( op.entity, result );
	    axDeselect();
	  }
	  else axLinkToPtrArray( result, NULL );
	  cp++;
	}
      else if( ! nEntity ) break;

#ifdef DEBUG
      axCurrent();
#endif

    }
  
  if( nEntity ) axPathValue( op.entity, result );
  else axLinkToPtrArray( result, NULL );

  return nEntity;
}

char *convertToAxOperation( char *logicalPath, struct axOperation *axOp )
{
    char buf[ AXENTITYLEN ];
    char *cp, *cp2;
    int n;

    cp = logicalPath;

    /* reset operation details */
    memset( axOp->entity, '\0', sizeof( char ) * AXENTITYLEN );
    memset( axOp->refinement, '\0', sizeof( char ) * AXENTITYLEN );
    axOp->selection = 0;

    /* Get name of concept or property */
    n = 0;
    while( isalpha( *cp ) && n < AXENTITYLEN - 1 ) axOp->entity[ n++ ] = *cp++;
    axOp->entity[ n ] = '\0';

    /* See if we have to select a particular no or do a refinement */
    while( *cp == '[' )
    {
        cp++;
        while( isspace(*cp) ) cp++;

	if( *cp == '*' ) 
	  {
	    
	    if( axOp->selection )
	      {
                fprintf( stderr, "Already set selection no!\n" );
                return NULL;
	      }
            axOp->selection = -2;
            cp++;
	  }
	
	/* get last in list */
	else if( *cp == '$' )
	  {
            if( axOp->selection )
	      {
                fprintf( stderr, "Already set selection no!\n" );
                return NULL;
	      }
            axOp->selection = -1;
            cp++;
	  }
	
        /* get specified no in list */
        else if( isdigit( *cp ) )
        {
            if( axOp->selection )
            {
                fprintf( stderr, "Already set selection no!\n" );
                return NULL;
            }

            n = 0;
            while( isdigit( *cp ) && n < AXENTITYLEN ) buf[ n++ ]= *cp++;
            buf[ n ] = '\0';

            axOp->selection = atoi( buf );

            if( axOp->selection < 1 )
            {
                fprintf( stderr, "Invalid selection no %d\n", axOp->selection );
                return NULL;
            }
        }

        /* get refinement to do */
        else if( isalpha( *cp ) )
        {
            if( *(axOp->refinement) )
            {
                fprintf( stderr, "Already set a refinement!\n" );
                return NULL;
            }

            n = 0;
            while( *cp && *cp != ']' && n < AXENTITYLEN - 1 )
            {
                if( *cp != '\'' ) axOp->refinement[ n++ ] = *cp++;
                else            
                {
                    /* AgentX doesn't seem to handle this properly */
                    /* axOp->refinement[n++]='\"'; */
                    cp++;
                }
            }
            axOp->refinement[ n ] = '\0';
        }

        else
        {
            fprintf( stderr, "Unknown refinement or selection!\n" );
            return NULL;
        }

        while( isspace(*cp) ) cp++;

        /* Check ] is present */
        if( *cp != ']' )
        {
            fprintf( stderr, "Malformed value selector!\n" );
            return NULL;
        }

        /* See if there is an additional selection/refinement */
        cp2 = cp + 1;
        while( isspace( *cp2 ) ) cp2++;
        if( *cp2 == '[' ) cp = cp2;

        /* don't increment pointer here so that we can distinguish
         * between error and reaching end of logical path. I know, I know... */
    }

#ifdef DEBUG
 printf( "convertToAxOperation: entity=%s\n", axOp->entity );
 printf( "convertToAxOperation: refinement=%s\n", axOp->refinement );
 printf( "convertToAxOperation: selection=%d\n", axOp->selection );
#endif

    return cp;
}

int doAxOperation( struct axOperation *axOp, char *cp, struct _axArray **result )
{
    int i , nEntity, cdepth, ndepth;

#ifdef DEBUG
    printf( "doAxOperation: entity: %s\n", axOp->entity );
    if( *( axOp->refinement ) ) printf( "doAxOperation: refinement: %s\n", axOp->refinement );
    if( axOp->selection ) printf( "doAxOperation: selection: %d\n", axOp->selection );
#endif

    /* Select entity */
    if( ( nEntity = axSelect( axOp->entity ) ) < 0 )
      fprintf( stderr, "Error selecting: %s\n", axOp->entity );
    if( nEntity <= 0 ) return nEntity;


    /* Process refinement string */
    if( *(axOp->refinement) )
    {
      if( ( nEntity = axRefine( axOp->refinement ) ) < 0 )
	fprintf( stderr, "Error refining: %s\n", axOp->refinement );
      if( nEntity <=0 ) return nEntity;
    }

    /* handle wildcards */

    /* Select specific list member */
    if( axOp->selection )
    {
      if( axOp->selection == -2 ){
	
	/* find the current depth */
	cdepth = axDepth();
	
	for( i = 1; i <= nEntity; i++ ){
	  
	  axPathNav( cp, result );
	  ndepth = axDepth();

	  while( (ndepth--) > cdepth ) axDeselect();
	  axSelectNext();
	  
	}
	
      }

      else{
        /* get last entity in list */
        if( axOp->selection == -1 ) axOp->selection = nEntity;
	
        else if( axOp->selection > nEntity )
	  {
            fprintf( stderr, "Invalid selection no: %d > %d (no of entities)\n", axOp->selection, nEntity );
            return -1;
	  }
	
        /* Move to specified entity */
        if( axSelectNo( axOp->selection ) < 0 )
	  {
            fprintf( stderr, "Error selecting entity %s no %d!\n", axOp->entity, axOp->selection );
            return -1;
	  }
      }
    }
    
    return nEntity;
}

/**************************************************************************/

