//***********************************************************************
//
//	Name:			CIMatFull.h
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




#ifndef __CIMatFull_h
#define __CIMatFull_h

#include "../../../../config.h"

#include "CIMat.h"




class CIMatFull : public CIMat {
public:
	CIMatFull(INT refDim, INT dim);
	CIMatFull(const MatrixType *p, INT refDim, INT dim);

	CIMatFull(istream &);
	~CIMatFull();

	CIMatFull(const CIMatFull &);
	CIMatFull(const CIMatFull &, INT dim);
	CIMatFull(const CIMatFull &, const MatSel &);
	CIMatFull & operator = (const CIMatFull &);


	MatrixType * getP();

	MatrixType &	operator()(INT i, INT j);
	MatrixType	operator()(INT i, INT j) const;

//----------------------------------------------------------------------------
//	virtual functions from DavidsonMatrix



	// get total matrix (attention: probably VERY large)
	void	getTotalMatrix(MatrixType *) const;
	void	getTotalMatrix(MatrixStorage<MatrixType, VectorType> &) const;


protected:
	INT	ind(INT i, INT j) const;

MatrixType	*p;
};





inline
MatrixType & CIMatFull::operator()(INT i, INT j)
{	return p[ind(i, j)];	}

inline
MatrixType	CIMatFull::operator()(INT i, INT j) const
{	return p[ind(i, j)];	}



inline
MatrixType * CIMatFull::getP()
{	return	p;	}


inline
INT	CIMatFull::ind(INT i, INT j) const
{	if ( i>j )
		return	i*(i+1)/2+j;
	else
		return	j*(j+1)/2+i;
}


#endif
