//***********************************************************************
//
//	Name:			InternalCoord.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			07. Feb 1998
//
//***********************************************************************




#ifndef __InternalCoord_h
#define __InternalCoord_h

#include "../../../config.h"


#include "Molecule.h"

#include "../../Container/String.h"

//FD class	ostream;
//FD class	istream;
class	GrepAwk;


class	InternalCoord {
public:
	InternalCoord();
	InternalCoord(const Molecule &);
	InternalCoord(istream &);
	InternalCoord(GrepAwk &);
	~InternalCoord();
	
	operator Molecule();
	
struct TInternal {
		String	Label;
		double	r;
		double	phi;
		double	theta;
		INT		rOpt;
		INT		phiOpt;
		INT		thetaOpt;
		INT		rel0;
		INT		rel1;
		INT		rel2;
	};
	
	void	writeToStream(ostream &) const;

	friend ostream & operator << (ostream &, const InternalCoord &);


private:
	void	read(GrepAwk &ga);

	Point	getCartCoord(Point p0, Point p1, Point p2, double r, double phi, double theta) const;

TInternal	*internal;
INT	n;
};







#endif
