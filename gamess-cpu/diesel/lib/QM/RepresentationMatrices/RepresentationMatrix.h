//***********************************************************************
//
//	Name:			RepresentationMatrix.h
//
//	Description:			
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __REPRESENTATIONMATRIX_H
#define __REPRESENTATIONMATRIX_H

#include "../../../config.h"

#include <iostream>

#include "CIVectors.h"
#include "../IntegralContainer/IntegralType.h"

template <class MatrixType> class RepresentationMatrix;
template <class MatrixType> ostream & operator << (ostream & s, const RepresentationMatrix<MatrixType> &);

template <class MatrixType>
class RepresentationMatrix {
public:
	
enum	MatrixDense {
		unity,
		diagonal,
		tridiagonal,
		sparse,
		full,
		sameSparse,
		lowerTriangle
	};
	
//------------------------------------------------------------------------

	// initialized static calling tables
	RepresentationMatrix();
	RepresentationMatrix(MatrixType *, INT r, INT c, INT simplify = 1);
	RepresentationMatrix(INT n, MatrixType *, INT r, INT c, MatrixDense dense);
//	RepresentationMatrix(MatrixType **, INT n, INT r, INT c);
//	RepresentationMatrix(RepresentationMatrix<MatrixType> &, MatrixType *);
	~RepresentationMatrix();


	// copy constructor
	RepresentationMatrix(const RepresentationMatrix<MatrixType> &);

//------------------------------------------------------------------------
struct TMPos {					// pointer to position of 
	INT	r, c;					// sparse matrix elements
	};

struct TMatDim {
	INT	rows;					// number of rows of matrices
	INT	cols;					// number of columns of matrices
	};
	
//----------------------------------------------------------------------	
	
	INT	getNSparse() const;
	const TMPos *getMPosP() const;
	
	INT	getNumberOfRows() const;
	INT	getNumberOfColumns() const;

	MatrixDense	getDense() const;
	void	setDense(MatrixDense);

	const MatrixType	*getP() const;
	
	INT	getOccupiedMemory() const;
	INT	getNElements() const;
	INT	getNMats() const;


		
//------------------------------------------------------------------------

	friend ostream & operator << <MatrixType> (ostream & s, const RepresentationMatrix<MatrixType> &);	

//------------------------------------------------------------------------

private:
	void	round(MatrixType *p, INT n);



//------------------------------------------------------------------------

MatrixDense	dense;				// type of matrix

TMatDim	matDim;					// dimension of matrix

INT	nSparse;					// number of entries in sparse matrices

TMPos	*mpos;

MatrixType	*elements;			// pointer to matrix elements
INT	nElements;					// number of individual elements
INT	nMats;						// number of matrices (for P=5: number!=1)
};

template <class MatrixType>	ostream & operator << (ostream & s, const RepresentationMatrix<MatrixType> &);	



template <class MatrixType>
inline
INT	RepresentationMatrix<MatrixType>::getNSparse() const
{	return nSparse;	}

template <class MatrixType>
inline
INT	RepresentationMatrix<MatrixType>::getNElements() const
{	return nElements;	}

template <class MatrixType>
inline
INT	RepresentationMatrix<MatrixType>::getNMats() const
{	return nMats;	}

template <class MatrixType>
inline
const typename RepresentationMatrix<MatrixType>::TMPos *RepresentationMatrix<MatrixType>::getMPosP() const
{	return mpos;	}

template <class MatrixType>
inline
INT	RepresentationMatrix<MatrixType>::getNumberOfRows() const
{	return matDim.rows;	}

template <class MatrixType>
inline
INT	RepresentationMatrix<MatrixType>::getNumberOfColumns() const
{	return matDim.cols;	}

template <class MatrixType>
inline
typename RepresentationMatrix<MatrixType>::MatrixDense	RepresentationMatrix<MatrixType>::getDense() const
{	return dense;	}

template <class MatrixType>
inline
void	RepresentationMatrix<MatrixType>::setDense(MatrixDense d)
{	dense = d;	}

template <class MatrixType>
inline
const MatrixType	*RepresentationMatrix<MatrixType>::getP() const
{	return	elements;	}



#endif
