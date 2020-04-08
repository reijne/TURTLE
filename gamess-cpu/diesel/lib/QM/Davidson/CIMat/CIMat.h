//***********************************************************************
//
//	Name:			CIMat.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			17. Dec 1998
//
//***********************************************************************

#ifndef __CIMat_h
#define __CIMat_h

#include "../../../../config.h"

#include "../DavidsonMatrix.h"

#include <stdlib.h>
#include <iostream>
using std::ostream;
using std::istream;

/* FD class ostream; */
/* FD class istream; */
class MatSel;

template <class MatrixType, class VectorType> class MatrixStorage;

typedef double MatrixType;
typedef double VectorType;

class CIMat : public DavidsonMatrix<MatrixType, VectorType> {
public:
	CIMat(INT refDim, INT dim);
	CIMat(const CIMat &mat, INT totalDim);
    	CIMat(istream &);
	~CIMat();

	CIMat(const CIMat &mat);
	CIMat & operator = (const CIMat & mat);

	void	setRandomCI(double dens,
		double diag, double diagDispl, double nonDiagDispl, INT randSeed);

	virtual MatrixType &	operator()(INT i, INT j) = 0;
	virtual MatrixType	operator()(INT i, INT j) const = 0;


//----------------------------------------------------------------------------
//	virtual functions from DavidsonMatrix

	// get diagonal of matrix
	void	getDiagonal(MatrixType *) const;
	

	// get reference part of matrix
	void	getRefMatrix(MatrixType *) const;

	// get index of nth reference in total matrix
	INT	getRefIndex(INT i) const;


	
//----------------------------------------------------------------------------

	void	chooseRefMat();


	friend ostream & operator << (ostream &, const CIMat &);

protected:
	void	mult(const double *x, double *y, INT n) const;


MatrixType *refMat;	
INT	*index;
};




#endif
