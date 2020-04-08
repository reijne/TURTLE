
#include<axArray.h>

#define AXELELEN                    750     /* max length of AxTransform elements */
#define AXENTITYLEN                 100     /* max length of AgentX concepts/properties */
#define ERRLEN                      200     /* max length of error string */

#ifndef ERROR
#define ERROR {retval=-1;goto cleanup;}
#endif

struct axOperation
{
  char entity[AXENTITYLEN];
  char refinement[AXENTITYLEN];
  int selection;
};

int axPath( char *path, char ***strArray );
int axPathNav( char *path, struct _axArray **result );
char *convertToAxOperation( char *logicalPath, struct axOperation *axOp );
int doAxOperation( struct axOperation *axOp, char *cp, struct _axArray **result );
