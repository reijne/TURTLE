//***********************************************************************
//
//	Name:			CIMatSparse.h
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

#ifndef __CIMatSparse_h
#define __CIMatSparse_h

#include "CIMat.h"
#include "../../MRCIMatrix/MatrixStorage.h"

#include <iostream>

class CIMatSparse : 
	public CIMat,
	protected MatrixStorage<MatrixType, VectorType> {
public:
	CIMatSparse(INT refDim, INT dim);
	CIMatSparse(istream &);
	~CIMatSparse();

	CIMatSparse(const CIMatSparse &);
	CIMatSparse(const CIMatSparse &, INT dim);
	CIMatSparse(const CIMatSparse &, const MatSel &);
	CIMatSparse & operator = (const CIMatSparse &);



	void	append(INT i, INT j, MatrixType v);
	MatrixType &	operator()(INT i, INT j);
	MatrixType	operator()(INT i, INT j) const;


//----------------------------------------------------------------------------
//	virtual functions from DavidsonMatrix

	// get total matrix (attention: probably VERY large)
	void	getTotalMatrix(MatrixType *) const;
	void	getTotalMatrix(MatrixStorage<MatrixType, VectorType> &) const;
	
//----------------------------------------------------------------------------

	void	getDiagonal(MatrixType *) const;


protected:
	INT	ind(INT i, INT j) const;
	void	mult(const double *x, double *y, INT n) const;

private:
static MatrixType zero;	
};




inline
void	CIMatSparse::append(INT i, INT j, MatrixType v)
{
	addEntry(i, j, v);
}


#endif
