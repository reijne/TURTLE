
/* SWIG file for the AgentX Perl wrapper */
 
%include<cstring.i>

%inline %{
#include<agentx.h>
%}

%cstring_output_allocate(char **, );						     

%typemap(in,numinputs=0) char *** ( char **temp ){

  temp = NULL;
  $1 = &temp;

}

%typemap(argout) char ***{

  int i = 0, sizeinc;
  SV **maxptr;
  
  maxptr = (*Perl_Tstack_max_ptr(((PerlInterpreter *)pthread_getspecific((*Perl_Gthr_key_ptr(((void *)0)))))));
  sizeinc = ( result - ( maxptr - sp ) );
  sp = Perl_stack_grow(((PerlInterpreter *)pthread_getspecific((*Perl_Gthr_key_ptr(((void *)0))))), sp,sp,sizeinc);
  
  if( temp$argnum ) for( ; i < result; i++ ){
    ST( i + argvi ) = SWIG_FromCharPtr( temp$argnum[ i ] );
  }
  
  argvi += i;
}

%include<agentx.h>
