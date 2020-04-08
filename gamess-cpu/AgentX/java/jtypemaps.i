
/* SWIG file for the AgentX Java wrapper */

%include<various.i>
%apply char **STRING_OUT { char ** };

%typemap(jni) char ***STRING_ARRAY "jobjectArray"

%typemap(in) char ***STRING_ARRAY ( char **temp ) {
  
  if ( ! $input ) {
    SWIG_JavaThrowException( jenv, SWIG_JavaNullPointerException, "array null" );
    return $null;
  }
  if ( JCALL1( GetArrayLength, jenv, $input ) == 0 ) {
    SWIG_JavaThrowException(jenv, SWIG_JavaIndexOutOfBoundsException, "Array must contain at least 1 element");
    return $null;
  }
  
  temp = NULL;
  $1 = &temp;
  
}

%typemap(argout) char ***STRING_ARRAY {
  
  int i, len = 0;
  jstring tmpStr = NULL;
  jclass strClass;
  jarray strArray;
  
  if( *$1 ) while ( (*$1)[ len++ ] );
  
  strClass = JCALL1( FindClass, jenv, "java/lang/String" );
  strArray = JCALL3( NewObjectArray, jenv, len, strClass, NULL );
  
  for( i = 0; i < len; i++ ) {
    tmpStr = JCALL1( NewStringUTF, jenv, (*$1)[ i ] );
    JCALL3( SetObjectArrayElement, jenv, strArray, i, tmpStr );
    JCALL1( DeleteLocalRef, jenv, tmpStr );
  }
  
  JCALL3( SetObjectArrayElement, jenv, $input, 1, strArray );
  
}

%apply char ***STRING_ARRAY { char *** };

%include<agentx.h>

