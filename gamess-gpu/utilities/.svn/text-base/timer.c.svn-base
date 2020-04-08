#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>    /* for gettimeofday */

long long walltime();

main()
{
  printf("%lld\n",walltime());
  exit(0);
}

long long walltime()
{
  static struct timeval tp;
  long long elp;
#if defined _SYSTYPE_SVR4
  int tzp;
#else
#if ! defined NEC
  struct timezone tzp;
#endif
#endif
   int t;
#if defined NEC
    gettimeofday(&tp);
#else
    gettimeofday(&tp,&tzp);
#endif
    elp =  1*(tp.tv_sec-1111770132) + (int) ( tp.tv_usec / 1000000.0 );
    return elp;
}
