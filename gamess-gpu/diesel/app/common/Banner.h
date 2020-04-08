//***********************************************************************
//
//	Name:			Banner.h
//
//	Description:	produces banner in output
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.08.1997
//
//***********************************************************************

#ifndef __Banner_h
#define __Banner_h

#include "../../config.h"

#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <string.h>


void	center(const char *s, INT w)
{
INT     w3 = strlen(s);
INT	w1 = (w-2-w3)/2;
INT	w2 = (w-2-w3-1)/2;
	if ( w1<0 )
		w1 = 1;
	if ( w2<0 )
		w2 = 1;
	cout << "*" << std::setw(w1) << " " << s << std::setw(w2) << " " << "*" << endl;
}


void	MakeBannerTop(INT w, const char *title)
{
	cout << 
"*******************************************************************************" << endl;
	center("", w);
	center(title, w);
	center("", w);
	center("part of DIESEL-MR-CI", w);
	center("", w);
}



void	MakeBannerBottom(INT w)
{
	center("", w);
	center("", w);
	center("by", w);
	center("Michael Hanrath and Bernd Engels", w);
	center("", w);
	center("SGA calculation: Volker Pleß", w);
	center("RI integral evaluation: Jan Franz", w);
	center("", w);
	center("Institute of Theoretical Chemistry", w);
	center("University of Bonn", w);
	center("Germany", w);
	center("", w);
	center("", w);
	center("", w);
	center("please cite as", w);
	center("M. Hanrath, B. Engels:", w);
	center("\"New algorithms for an individually selecting MR-CI program\",", w);
	center("Chemical Physics 225 (1997) 197-202", w);
	center("", w);
	center("", w);
	center("", w);
	cout << 
"*-----------------------------------------------------------------------------*" << endl;
	center("", w);
	center("", w);
	CONFIGURED_BANNER
	COMPILED_BANNER
	center("", w);
	center("", w);
	cout << 
"*******************************************************************************" << endl;
	cout << endl << endl << endl;	

char	buf[1000];
time_t s_time = time(NULL);

	cout << 
"*******************************************************************************" << endl;
	cout << endl;
	cout << "\t\t\tuser      : " << getenv("LOGNAME") << endl;
	gethostname(buf, 1000);
	cout << "\t\t\thostname  : " << buf << endl;
	cout << "\t\t\ttime/date : " << ctime(&s_time) << endl;
	cout << 
"*******************************************************************************" << endl;
	cout << endl << endl << endl;	
}


#endif

