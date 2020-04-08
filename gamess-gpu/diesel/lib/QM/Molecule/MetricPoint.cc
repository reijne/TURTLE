//***********************************************************************
//
//	Name:			MetricPoint.cc
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




#include "MetricPoint.h"

#define EPS 1E-6


MetricPoint::MetricPoint(INT nPhi, INT potPhi)
		: Point()
{
	setX(cos(2*M_PI*potPhi/nPhi));
	setY(sin(2*M_PI*potPhi/nPhi));
	setZ(0);
}

		



MetricPoint::MetricPoint(INT nPhi, INT potPhi, double Theta)
	: Point()
{
	setMetric(nPhi, potPhi, Theta);
}


void	MetricPoint::setMetric(INT nPhi, INT potPhi, double Theta)
{
double	h = sin(Theta);

	if ( nPhi==0 )
	{	setX(h);
		setY(0);
	}
	else
	{	setX(h*cos(2*M_PI*potPhi/nPhi));
		setY(h*sin(2*M_PI*potPhi/nPhi));
	}
	setZ(sqrt(1-h*h));
}

		
INT operator==(const MetricPoint & p1, const MetricPoint & p2)
{
INT	n1, n2, pot1, pot2;
double	theta1, theta2;
	p1.getType(theta1, n1, pot1);
	p2.getType(theta2, n2, pot2);

	return  (fabs(theta1-theta2)<EPS) && 
		(n1==n2) && (pot1==pot2);
}

		
INT operator<(const MetricPoint & p1, const MetricPoint & p2)
{
INT	n1, n2, pot1, pot2;
double	theta1, theta2;
	p1.getType(theta1, n1, pot1);
	p2.getType(theta2, n2, pot2);

INT	gle=p1.calcGle(theta1, theta2);

	if ( gle<0 )
		return 1;
		
	if ( gle>0 )
		return 0;
		
	if ( n1<n2 )
		return 1;

	if ( n1>n2 )
		return 0;
		
	return pot1<pot2;
}


INT operator<=(const MetricPoint & p1, const MetricPoint & p2)
{
INT	n1, n2, pot1, pot2;
double	theta1, theta2;
	p1.getType(theta1, n1, pot1);
	p2.getType(theta2, n2, pot2);

INT	gle=p1.calcGle(theta1, theta2);

	if ( gle<0 )
		return 1;
		
	if ( gle>0 )
		return 0;
		
	if ( n1<n2 )
		return 1;

	if ( n1>n2 )
		return 0;
		
	return pot1<=pot2;
}


INT	MetricPoint::cmpNormal(MetricPoint p2) const
{
INT	n[2], pot[2];
double	theta[2];
	getType(theta[0], n[0], pot[0]);
	p2.getType(theta[1], n[1], pot[1]);

	

	for ( INT i=0 ; i<2 ; i++ )
	{	if ( theta[i]>M_PI/2 ) 	// theta > pi/2
		{	theta[i] = M_PI - theta[i];
			pot[i] += n[i]/2;
			pot[i] %= n[i];
		}
		else
		if ( fabs(theta[i]-M_PI/2)< EPS )		// theta = pi/2
			if ( n[i]/2 )
				if ( pot[i]>n[i]/2 )
					pot[i] -= n[i]/2;

		if ( !pot[i] )
			n[i] = 0;
	}

INT	gle=calcGle(theta[0], theta[1]);

	if ( gle<0 )
		return -1;
		
	if ( gle>0 )
		return 1;

	if ( n[0]<n[1] )
		return -1;

	if ( n[0]>n[1] )
		return 1;
		
	if ( pot[0]<pot[1] )
		return -1;

	if ( pot[0]>pot[1] )
		return 1;
		
	return 0;
}


MetricPoint	MetricPoint::operator-()
{
	return MetricPoint(-getX(), -getY(), -getZ());
}

void	MetricPoint::getType(
				double & Theta,
				INT & nPhi, INT & potPhi
				) const
{
double	x = getX();
double	y = getY();
double	z = getZ();

double	r = sqrt(x*x + y*y + z*z);
double	phi;


	if ( r==0 )
		Theta = phi = 0;
	else
	{	//cout << "z=" << z << ", r=" << r << endl;
		Theta = acos(z/r);
		phi = atan2(y, x)/(2*M_PI);
		if ( Theta<0 )
			Theta += 1;

		if ( phi<0 )
			phi += 1;
	}
	potPhi = 0;
	nPhi = 1;

	if ( fabs(Theta) < EPS )
	{	Theta = 0;
		return;
	}	
		
	getSmallestRational(phi, nPhi, potPhi, EPS);
	if ( !potPhi )
		nPhi = 0;
	if ( nPhi == potPhi )
		nPhi = potPhi = 0;
//	cout << *this;
//	cout << Theta << endl;
//	cout << phi << " " << nPhi << " " << potPhi << endl;

}


INT	test()
{
  return 0;
}

/*
void	MetricPoint::getType(
				INT & nTheta, INT & nPhi,
				INT & potTheta, INT & potPhi,
				INT & isDiag,
				INT N
				) const
{
	getType(nTheta, nPhi, potTheta, potPhi);
	if ( (nPhi==2*N) && (potPhi & 1) )
		isDiag = 1;
	else
		isDiag = 0;
}
*/


void	MetricPoint::getSmallestRational(
	double alpha, INT & nalpha, INT & potalpha, double eps) const
{
INT	i;
	for ( i=1 ; i<MaxPot ; i++ )
		if ( fabs(alpha*i-floor(alpha*i+0.5))<eps )
			break;
	if ( i<MaxPot )
	{	nalpha = i;
		potalpha = (INT) floor(alpha*i+0.5);
	}
	else
	{	nalpha = 0;
		potalpha = 0;
	}	
}


INT	MetricPoint::calcGle(double theta1, double theta2) const
{
double	diff = theta1-theta2;

	if ( fabs(diff)<EPS )
		return 0;
		
	if ( fabs(theta1)<EPS )
		return -1;
		
	if ( fabs(theta2)<EPS )
		return 1;
		
	if ( fabs(theta1-M_PI/2)<EPS )
		return -1;
		
	if ( fabs(theta2-M_PI/2)<EPS )
		return 1;
		
	if ( diff<0 )
		return -1;
	else
		return 1;
}
