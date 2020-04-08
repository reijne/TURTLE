/*$Id: ieeetrap.c,v 1.1.1.1 2000-10-26 16:29:41 psh Exp $*/
#include <floatingpoint.h>
#include <stdio.h>
#include <signal.h>

static void catchit()
{
  printf("!!  Floating point interrupt caught  !!\n");
  fflush(stdout);
  (void) signal(SIGIOT, SIG_DFL);
  abort();
}

void ieeetrap_()
{
 (void) ieee_handler("set","inexact", SIGFPE_IGNORE);
 (void) ieee_handler("set","underflow", SIGFPE_IGNORE);
 (void) ieee_handler("set","invalid", SIGFPE_IGNORE);

}
