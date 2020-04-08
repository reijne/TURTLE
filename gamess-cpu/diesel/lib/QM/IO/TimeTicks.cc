//***********************************************************************
//
//	Name:			TimeTicks.cc
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


#include "TimeTicks.h"


#include <string.h>
#include <iomanip>


TimeTicks::TimeTicks()
{
	clear();
}

void	TimeTicks::clear()
{
	memset(&sTms, 0, sizeof(tms));
	memset(&eTms, 0, sizeof(tms));
	sWall = 0;
	eWall = 0;
}



void	TimeTicks::start()
{
	times(&sTms);
	time(&sWall);
}


void	TimeTicks::stop()
{
tms	_tms;
	times(&_tms);
	eTms.tms_utime = _tms.tms_utime - sTms.tms_utime;
	eTms.tms_stime = _tms.tms_stime - sTms.tms_stime;
	eTms.tms_cutime = _tms.tms_cutime - sTms.tms_cutime;
	eTms.tms_cstime = _tms.tms_cstime - sTms.tms_cstime;

time_t	h;
	time(&h);
	eWall = h - sWall;
}



INT	TimeTicks::getUserTicks() const
{	return eTms.tms_utime;	}

INT	TimeTicks::getSysTicks() const
{	return eTms.tms_stime;	}

INT	TimeTicks::getUserChildrenTicks() const
{	return eTms.tms_cutime;	}

INT	TimeTicks::getSysChildrenTicks() const
{	return eTms.tms_cstime;	}


double	TimeTicks::getUser() const
{	return ((double) eTms.tms_utime)/CLOCKS_PER_SEC;		}

double	TimeTicks::getSys() const
{	return ((double) eTms.tms_stime)/CLOCKS_PER_SEC;		}

double	TimeTicks::getUserChildren() const
{	return ((double) eTms.tms_cutime)/CLOCKS_PER_SEC;	}

double	TimeTicks::getSysChildren() const
{	return ((double) eTms.tms_cstime)/CLOCKS_PER_SEC;	}


time_t	TimeTicks::getWall() const
{	return eWall;	}




TimeTicks & TimeTicks::operator+=(const TimeTicks & tt)
{
	eTms.tms_utime += tt.getUserTicks();
	eTms.tms_stime += tt.getSysTicks();
	eTms.tms_cutime += tt.getUserChildrenTicks();
	eTms.tms_cstime += tt.getSysChildrenTicks();
	eWall += tt.getWall();
	return *this;
}


ostream& operator<<(ostream & s, const  TimeTicks & ticks)
{
INT	prec = s.precision();

	s << "(";
	s << std::setprecision(2) << "Wall: " << ticks.getWall() << "s, ";

	s << "user: " << ticks.getUser() << "s";
	if ( ticks.getUserChildren()>0 )
		s << "/" << ticks.getUserChildren() << "s";

	s << ", ";

	s << "sys: " << ticks.getSys() << "s";
	if ( ticks.getSysChildren()>0 )
		s << "/" << ticks.getSysChildren() << "s";
	s << ")";
	s.precision(prec);
	return s;
}
