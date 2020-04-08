
#define MAX_BUF_LEN  ( MAX_STR_LEN * 5 )
#define ERROR_SP     stderr

#ifndef ERROR
#define ERROR {retval=-1;goto cleanup;}
#endif

int axCheckBufLen( char *str, unsigned int size );
int axRegError( char *error );
int axClearErrors();
int axPrintErrors();
int axSetErrorStream( FILE *es );
int axGetError( char **error );
int axGetErrors( char ***errors );
