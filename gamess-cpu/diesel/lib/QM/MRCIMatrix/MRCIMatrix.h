//***********************************************************************
//
//	Name:			MRCIMatrix.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.06.1996
//					21.02.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MRCIMATRIX_H
#define __MRCIMATRIX_H

#include "../../../config.h"


#include <math.h>


#include "../Configuration/DiffConf.h"
#include "../Configuration/TableCase.h"

#include "../Davidson/DavidsonMatrix.h"

#include "CICalculation.h"

#include "../RepresentationMatrices/CIVectors.h"


class NExternalsDiag;
class TupelStructureDiag;

class DiskBuffer;

template <class MatrixType, class VectorType> class MatrixStorage;
template <class Key, class ContainedObject> class VarSizeReadOnlyCache;


template <class MatrixType, class VectorType>
class MRCIMatrix : 
	public CICalculation<MatrixType, VectorType>,
	public DavidsonMatrix<MatrixType, VectorType> {
public:

	//	for hamilton matrix multiplication
	MRCIMatrix(NExternalsDiag	*,
		Fort31File _f31,
		INT cacheEntries = 2048, INT cacheMemory = (1 << 20),
		INT numberOfSlaves = 0);

	//	for density matrix calculation
	MRCIMatrix(NExternalsDiag *, NExternalsDiag *,
		Fort31File _f31,
		INT cacheEntries = 2048, INT cacheMemory = (1 << 20),
		INT numberOfSlaves = 0);

	//	for reference density matrix calculation
	MRCIMatrix(const CICalculation<MatrixType, VectorType> &, NExternalsDiag *);

	~MRCIMatrix();



	const NExternalsDiag *
		getNExternalsDiag() const;



	

	INT	getRefIndex(INT i) const;


	//----------------------------------------------------------------	
	//	virtual functions from DavidsonMatrix

	//	perform y = A*x
	void	mult(DiskBuffer *x, DiskBuffer *y, INT start, INT end) const
	{	((MRCIMatrix *) this)->_mult(x, y, start, end); }
	
	//	get diagonal
	void	getDiagonal(MatrixType *p) const
	{	((MRCIMatrix *) this)->_getDiagonal(p); }
	
	//	get ci-matrix of reference space
	void	getRefMatrix(MatrixType *pRefMatrix) const;
	void	getReferenceDensity(CIVectors<VectorType> & x, VectorType *densMat) const;
	void	getDensity(CIVectors<VectorType> & x, VectorType *densMat) const;

	//	get complete ci-matrix (full matrix, v e r y large!)
	//	only for test purposes
	void	getTotalMatrix(MatrixType *pCIMatrix) const;
	
	void	getTotalMatrix(MatrixStorage<MatrixType, VectorType> &m) const;


	//----------------------------------------------------------------	

	LONG_LONG_INT estimateNonZeroElements(INT maxBreak) const;

	//----------------------------------------------------------------	

	//	calculate density matrix
	void	dens(const CIVectors<VectorType> * x, const CIVectors<VectorType> * y,
		VectorType **densMat) const;

	//----------------------------------------------------------------	

	INT	getNumberOfSlaves() const;

	CIVectors<VectorType> *	getXVector();
	CIVectors<VectorType> *	getYVector();

	//----------------------------------------------------------------	

	VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*getCache() const;

	//----------------------------------------------------------------	

private:
//============================================================================

struct	MultData {
enum CalcMode { _Multiplication, _Density, _MatrixStorage, _PointerStorage };

	MultData(CIVectors<VectorType> * _x, CIVectors<VectorType> * _y)
	{
		calcMode = _Multiplication;
		x = _x;
		y = _y;
		stored = NULL;
		densMat = NULL;
		HMat = NULL;
	}
	MultData(MatrixStorage<MatrixType, VectorType> * _stored)
	{
		calcMode = _MatrixStorage;
		x = NULL;
		y = NULL;
		stored = _stored;
		densMat = NULL;
		HMat = NULL;
	}
	MultData(const CIVectors<VectorType> * _x, const CIVectors<VectorType> * _y, 
		VectorType	**_densMat)
	{
		calcMode = _Density;
		x = (CIVectors<VectorType> *)_x;
		y = (CIVectors<VectorType> *)_y;
		stored = NULL;
		densMat = _densMat;
		HMat = NULL;
	}
	MultData(MatrixType	*_HMat)
	{
		calcMode = _PointerStorage;
		x = NULL;
		y = NULL;
		stored = NULL;
		densMat = NULL;
		HMat = _HMat;
	}
	CalcMode	calcMode;
	CIVectors<VectorType> *x;
	CIVectors<VectorType> *y;
	MatrixStorage<MatrixType, VectorType> *stored;	// stored hamilton matrix
	VectorType	**densMat;							// pointer to density matrices
	MatrixType	*HMat;								// pointer to hamilton matrix
	};

//============================================================================

	void	mult(const VectorType *x, VectorType *y, INT n) const {};

	void	_mult(DiskBuffer *x, DiskBuffer *y, INT start, INT end);
	void	_getDiagonal(MatrixType *p);
	void	_getRefMatrix(MultData data);
//	void	_getTotalMatrix(MatrixType *pCIMatrix);



	//	perform y = A*x
	void	calcHamilton(NExternalsDiag	*, NExternalsDiag *, MultData);
	
	//	perform y = A*x, parallel version
//	void	calcHamilton();
	
NExternalsDiag	*mrccA;				// configuration tree
NExternalsDiag	*mrccB;				// optional second configuration tree
									// (used for density matrix calculation)


INT	*RefSAFIndex;

VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*cache;

INT	numberOfSlaves;

};


template <class MatrixType, class VectorType>
inline
INT	MRCIMatrix<MatrixType, VectorType>::getRefIndex(INT i) const
{	return RefSAFIndex[i];	}

template <class MatrixType, class VectorType>
inline
INT	MRCIMatrix<MatrixType, VectorType>::getNumberOfSlaves() const
{	return	numberOfSlaves;	}


template <class MatrixType, class VectorType>
inline
VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	
	*MRCIMatrix<MatrixType, VectorType>::getCache() const
{	return cache;	}


#endif
