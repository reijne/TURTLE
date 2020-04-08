//***********************************************************************
//
//	Name:			StartUp.cc
//
//	Description:	install signal catcher and ensures a proper clean up
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.08.1997
//
//***********************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "../../config.h"

#include <iostream>

using namespace std;

#include "../../lib/QM/IO/TimeTicks.h"


extern TimeTicks globalTime;

void	cleanUp()
{
	globalTime.stop();
	cout << "total time: " << globalTime << endl;
	printf("cleaning up.\n");
}

void	catchSignals(INT i)
{
	globalTime.stop();
	cout << "total time: " << globalTime << endl;
#ifdef __HAS_STR_SIGNAL
	printf("catched signal %s\n", strsignal(i));
#elseif
	printf("catched signal %d\n", i);
#endif
	printf("cleaning up.\n");
}


void	StartUp()
{
	atexit(cleanUp);
/*	signal(SIGHUP, catchSignals);
	signal(SIGINT, catchSignals);
	signal(SIGQUIT, catchSignals);
	signal(SIGTERM, catchSignals);
	signal(SIGBUS, catchSignals);
	signal(SIGSEGV, catchSignals);
*/
}
