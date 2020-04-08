#include "../../../config.h"
//***********************************************************************
//
//	Name:			MOTrafo.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.05.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MOTrafo_h
#define __MOTrafo_h

#include "../../Math/MatrixVector/Matrix.h"

#include "../MO/MOIrReps.h"


/* FD class istream; */
/* FD class ostream; */

template <class T> class Matrix;
template <class T> class Vector;


class MOTrafo : public MOIrReps {
public:
	MOTrafo(const MOIrReps &, istream &s);
	MOTrafo(istream &s);
	~MOTrafo();
	
typedef double MOCoefType;

	operator Matrix<MOCoefType> () const;
	
	
	MOCoefType	operator () (INT sym, INT i, INT j) const;
	MOCoefType	operator () (INT i, INT j) const;


	void	transform(const Matrix<MOCoefType> &);
	void	setOccupationNumbers(const Vector<MOCoefType> &);
	Vector<MOCoefType>	getOccupationNumbers() const;

	ostream & writeToStream(ostream &) const;
	friend ostream & operator << (ostream &, const MOTrafo &);
	
private:
	MOTrafo(const MOTrafo &);
	MOTrafo & operator = (const MOTrafo &);


	void	unpack(MOCoefType *) const;
	void	pack(const MOCoefType *);
	
MOCoefType	**coefs;		//	array of pointers to coefficients
MOCoefType	**occNum;		//	array of pointers to occupation numbers

};



inline
MOTrafo::MOCoefType	MOTrafo::operator () (INT sym, INT i, INT j) const
{
	return coefs[sym][i*inIrRep[sym] + j];
}

inline
MOTrafo::MOCoefType	MOTrafo::operator () (INT i, INT j) const
{
	if ( getIrRep(i)!=getIrRep(j) )
		return 0;
IrRep	sym = getIrRep(i);
	return coefs[sym][(i-getStartMO(sym))*inIrRep[sym] + j - getStartMO(sym)];
}

#endif
