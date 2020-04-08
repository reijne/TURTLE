
#include<stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  /* load a supplementary document */
  int axGetUri( char* );
  
  /* load a data document */
  int axDataGetUri( char* );
  
  /* locate data sets relating to a concept */
  int axSelect( char* );
  
  /* select next data set */
  int axSelectNext();
  
  /* select previous data set */
  int axSelectPrev();
  
  /* drop the previously located data sets */
  int axDeselect();
  
  /* specify a base URI for query terms */
  int axBaseUri( char* );
  
  /* clean up parser */
  int axParserFinish();
  
  /* initialise parser */
  int axParserStart();
  
  /* refine data sets previously located */
  int axRefine( char* );
  
  /* set the data format
     1 = XML */
  int axFormat( int );
  
  /* display location history */
  int axCurrent();
  
  /* return a data element */
  char* axValue();
  
  /* set the size of the evaluation cache (in evaluations) */
  int axCache( int );
  
  /* select a specific data set */
  int axSelectNo( int );
  
  /* return a data element related to metadata */
  char* axAbout( char* );
  
  /* same as axSelect(); axValue(); axDeselect() */
  char* axSelectValue( char* );
  
  /* parse an AgentX control file*/
  int axControl( char *,char * );
  
  /* deprecated */
  int axSubstitute( char *, char ** );

  /* evaluate an axpath */
  int axPath( char *, char ***values );

  /* clear current error stack */
  int axClearErrors();

  /* write error stack to ERROR_SP (default = stderr) */
  int axPrintErrors();

  /* set ERROR_SP */
  int axSetErrorStream( FILE * );

  /* return the top level error */
  int axGetError( char **errorString );

  /* return the error stack */
  int axGetErrors( char ***errorStrings );

#ifdef __cplusplus
}
#endif
