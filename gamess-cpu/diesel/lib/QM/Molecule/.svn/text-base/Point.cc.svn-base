//***********************************************************************
//
//	Name:			Point.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03. Nov 1998
//
//***********************************************************************




#include "Point.h"


//template class Point;

Point::Point(double x, double y, double z) : Vector<double>(3)
{
	p[0] = x;
	p[1] = y;
	p[2] = z;
}

Point::Point(istream &s) : Vector<double>(3)
{
	s >> p[0];
	s >> p[1];
	s >> p[2];
}

Point::Point(const Point & p1) : Vector<double>(3)
{
	memcpy(p, p1.getP(), 3*sizeof(double));
}


Point & Point::operator=(const Point &p1)
{
	memcpy(p, p1.getP(), 3*sizeof(double));
	return *this;
}

Point & Point::operator=(const Vector<double> & p1)
{
	memcpy(p, p1.getP(), 3*sizeof(double));
	return *this;
}

Point::Point(const Vector<double> & v) : Vector<double>(3)
{
	memcpy(p, v.getP(), 3*sizeof(double));
}


/*
main()
{
Point	p = Point(1, 5, 3);
Vector<double>	v = Vector<double>(3, 1.0, 2.0, 3.0);

//	p.normalize2();
	v = v % p;

	cout << v << p;
}
*/


void	Point::writeToStream(ostream &s) const
{
	s << getX() << " " << getY() << " " << getZ();
}


ostream &	operator << (ostream &s, const Point & a)
{
	a.writeToStream(s);
	return s;
}
