

/* SWIG file for the AgentX Python wrapper */
 
%include<cstring.i>

%cstring_output_allocate(char **, );						     

%typemap(in,numinputs=0) char *** ( char **temp ){

  temp = NULL;
  $1 = &temp;

}

%typemap(argout) char ***{
  
  int i;
  char *str;

  if( temp$argnum ) for( i = 0; i < result; i++ ){
    %append_output( SWIG_FromCharPtr( temp$argnum[ i ] ) );
  }
}

%include<agentx.h>
