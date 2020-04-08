//***********************************************************************
//
//	Name:			MatrixStorage.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.08.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MatrixStorage_h
#define __MatrixStorage_h

//FDclass	ostream;

#include "../../../config.h"

#include "../IO/Fortran/FortranFileIO.h"


template <class MatrixType, class VectorType>
class MatrixStorage {
public:

	MatrixStorage(INT dim);
	~MatrixStorage();
	
	
	MatrixStorage(const MatrixStorage<MatrixType, VectorType> &);
	MatrixStorage<MatrixType, VectorType> & operator = (
		const MatrixStorage<MatrixType, VectorType> &);
	
	
	INT	getDim() const;
	INT	getNumberOfEntries() const;

	
	double	getSparsity() const;
	
	static LONG_LONG_INT getMemory(LONG_LONG_INT n);
	
	void	addEntry(INT r, INT c, MatrixType v);
	
	//	perform y = A*x
	void	mult(const VectorType *x, VectorType *y) const;
	void	mult(const VectorType *x, VectorType *y, INT n) const;
	

//	friend ostream & operator << <> (ostream &s, 
//		const MatrixStorage<MatrixType, VectorType> &);


protected:
	MatrixStorage();
	


INT	dim;						// dimension of matrix
INT	nEntries;					// number of entries	
struct TEntry {					// storageMode = NonZero
	MatrixType	v;					// value
	INT		row;				// row
	INT		col;				// column
} *p;
static const INT growthSizePot = 10;			// 
static const INT growthSize = 1 << 10;			// size of reallocation blocks
static const INT growthSizeMask = (1 << 10)-1;	// size mask of reallocation blocks


static const INT blocksizeInts = 1 << 12;
static const INT blocksizeBytes = sizeof(TEntry) * blocksizeInts;
TEntry	buffer[blocksizeInts];
INT	wo;
FortranFileIO	*fio;
};




#include <stdlib.h>
#include <math.h>

extern INT writeHamiltonOnly;

template <class MatrixType, class VectorType>
inline
INT	MatrixStorage<MatrixType, VectorType>::getDim() const
{	return	dim;	}


template <class MatrixType, class VectorType>
inline
INT	MatrixStorage<MatrixType, VectorType>::getNumberOfEntries() const
{	return	nEntries;	}


template <class MatrixType, class VectorType>
inline
void	MatrixStorage<MatrixType, VectorType>::addEntry(INT r, INT c, MatrixType v)
{
	if ( fabs(v)>1e-20 )
	{
		if ( !writeHamiltonOnly )
		{
			if ( !(nEntries & growthSizeMask) )
				p = (TEntry *) realloc(p, (nEntries+growthSize)*sizeof(TEntry));
			p[nEntries].row = r;
			p[nEntries].col = c;
			p[nEntries].v = v;
		}
		else
		{
			buffer[wo].row = r;
			buffer[wo].col = c;
			buffer[wo++].v = v;
			if ( wo>=blocksizeInts )
			{
				fio->write(buffer, blocksizeBytes);
				wo = 0;
			}
		}
		nEntries++;
	}
}


#endif
