//***********************************************************************
//
//	Name:			MetricPoint.h
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




#ifndef _METRICPOINT_H
#define _METRICPOINT_H

#include "../../../config.h"

#include "Point.h"


class MetricPoint : public Point {

public:
	MetricPoint(double x=0, double y=0, double z=1) :
		Point(x, y, z)
	{}
	
	MetricPoint(INT nPhi, INT potPhi);
//	MetricPoint(INT nPhi, INT potPhi, INT Diag);
	MetricPoint(INT nPhi, INT potPhi, double Theta);
//	MetricPoint(INT nPhi, INT potPhi, INT nTheta, INT potTheta, INT Diag);
		
	friend INT operator==(const MetricPoint &, const MetricPoint &);
	friend INT operator!=(const MetricPoint & s1, const MetricPoint & s2)
	{	return !(s1==s2);	}
	friend INT operator<(const MetricPoint &, const MetricPoint &);
	friend INT operator<=(const MetricPoint &, const MetricPoint &);
	friend INT operator>(const MetricPoint & s1, const MetricPoint & s2)
	{	return !(s1<=s2);	}
	
	friend INT operator>=(const MetricPoint & s1, const MetricPoint & s2)
	{	return !(s1<s2);	}

	MetricPoint operator-();

	INT	cmpNormal(MetricPoint p2) const;

	void	setMetric(INT nPhi, INT potPhi, double Theta);
	
	void	getType(
				double & Theta, 
				INT & nPhi, INT & potPhi
				) const;
				
	INT calcGle(double theta1, double theta2) const;
				
static const INT	MaxPot = 100;
private:
	void getSmallestRational(double alpha, INT & nalpha, INT & potalpha, 
		double eps) const;
};

#endif
