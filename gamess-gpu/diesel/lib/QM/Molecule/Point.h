//***********************************************************************
//
//	Name:			Point.h
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

#ifndef _POINT_H
#define _POINT_H

#include "../../../config.h"

#include "../../Math/MatrixVector/Vector.h"

#include <iostream>
using std::istream;
using std::ostream;

class Point : public Vector<double> {
public:

	Point(double x=0, double y=0, double z=0);
	Point(istream &s);

	Point(const Point &);

	Point & operator=(const Point &);
	Point & operator=(const Vector<double> &);
	Point(const Vector<double> &);



	
	double &	getX()
	{	return p[0];	}

	double &	getY()
	{	return p[1];	}

	double &	getZ()
	{	return p[2];	}
	
	double 	getX() const
	{	return p[0];	}

	double 	getY() const
	{	return p[1];	}

	double 	getZ() const
	{	return p[2];	}
	
	void	setX(double x)
	{	p[0] = x;	}
	
	void	setY(double y)
	{	p[1] = y;	}
	
	void	setZ(double z)
	{	p[2] = z;	}


	void	getCoordinates(double & x, double & y, double & z) const
	{	x = p[0];
		y = p[1];
		z = p[2];
	}
	
	void	setCoordinates(double & x, double & y, double & z)
	{	p[0] = x;
		p[1] = y;
		p[2] = z;
	}
	
	void	writeToStream(ostream &s) const;

	friend ostream &	operator << (ostream &s, const Point &);
	
protected:

};

#endif
