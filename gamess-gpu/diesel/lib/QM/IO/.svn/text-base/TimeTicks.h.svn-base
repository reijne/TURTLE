//***********************************************************************
//
//	Name:			TimeTicks.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27. Jan 1998
//
//***********************************************************************

#ifndef __TIMETICKS_H
#define __TIMETICKS_H

#include "../../../config.h"

#include <iostream>
using std::ostream;

#include <time.h>
#include <sys/times.h>

class TimeTicks {
public:
	TimeTicks();
	~TimeTicks() {}

//------------------------------------------------------------------

	void	start();
	void	stop();

	void	clear();

	INT	getUserTicks() const;
	INT	getSysTicks() const;
	INT	getUserChildrenTicks() const;
	INT	getSysChildrenTicks() const;

	double	getUser() const;
	double	getSys() const;
	double	getUserChildren() const;
	double	getSysChildren() const;

	time_t	getWall() const;

	TimeTicks & operator+=(const TimeTicks &);


//------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const  TimeTicks &);


private:
tms	sTms;
tms	eTms;
time_t	sWall;
time_t	eWall;
};




#endif
