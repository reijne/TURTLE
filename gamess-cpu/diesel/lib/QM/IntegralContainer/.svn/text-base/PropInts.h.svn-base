#include "../../../config.h"
//***********************************************************************
//
//	Name:			PropInts.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.05.1998
//
//
//
//
//
//***********************************************************************


#ifndef __PropInts_h
#define __PropInts_h

#include "../MO/MOIrReps.h"

// class istream;
// class ostream;

class MOTrafo;
class DensityMatrix;


class PropInts : public MOIrReps {
public:
typedef double Type;
enum	Operator { Mltpl1, Kinetic, OneHam, AngMom };
static char	*OperatorNames[];
	PropInts(Operator op, INT component, const char *filename = "ONEINT");
	~PropInts();


	void	transform(const MOTrafo &);

	Type	multDensity(const DensityMatrix &, bool anti) const;

	Type *	expand() const;
	void	collect(const Type *p) const;
	
	friend ostream & operator << (ostream &, const PropInts &);
	
private:

	INT	ind(INT i, INT j) const;
	INT	ind(IrRep irrep, INT i, INT j) const;

INT	nInts;					// number of integrals
INT	SymLbl;					// compressed product table
Type	*p;					// pointer to integrals
};



inline
INT	PropInts::ind(INT i, INT j) const
{
	if ( i>j )
		return (i+1)*i/2 + j;
	else
		return (j+1)*j/2 + i;	
}

inline
INT	PropInts::ind(IrRep irrep, INT i, INT j) const
{
	return i*inIrRep[irrep] + j;
}


#endif
