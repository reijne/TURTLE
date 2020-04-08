//***********************************************************************
//
//	Name:			MRMPH0MatElements.h
//
//	Description:	calculates Multi Reference Moller-Plesset
//					Hamilton Matrix elements
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.09.1998
//
//
//
//***********************************************************************

#ifndef __MRMPH0MatElements_h
#define __MRMPH0MatElements_h

#include "../../../config.h"

#include "RepresentationMatrices.h"


template <class DensMatType> class MRFockMatrix;



template <class MatrixType, class VectorType>
class MRMPH0MatElements : public RepresentationMatrices<MatrixType> {
public:
	// initialize static calling tables

	//	for hamilton matrix multiplication
	MRMPH0MatElements(
		const MRFockMatrix<VectorType> *,
		TwoIndexIntegralContainer *,
		FourIndexIntegralContainer *,
		double core,
		MatrixType E0);

	MRMPH0MatElements(TableKey key) :
		RepresentationMatrices<MatrixType>(key) {}
	~MRMPH0MatElements();

friend class RepDiag<MatrixType, VectorType>;	

//---------------------------------------------------------------------

	void	getMatrix(
		MatrixType *p,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const;

	void	getCaseP3Matrix(
		MatrixType *p,
		const DiffConf<MOType> & diffConf,
		INT	dK, INT R) const;
	
	void	getCaseP3Matrix(
		MatrixType *p,
		MOType a, MOType b) const;


	void	getCase14Matrix(
		MatrixType *p,
		const Configuration<MOType> & same) const;

	

//---------------------------------------------------------------------	

	void	multMatrix(
		VectorType *y, const VectorType *x,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const;

	void	multCaseP3Matrix(
		VectorType *y, const VectorType *x,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		INT	dK, INT R) const;
	
	void	multCaseP3Matrix(
		VectorType *y, const VectorType *x,
		INT SAFStartA, INT SAFStartB,
		MOType a, MOType b) const;


	void	multCase14Matrix(
		VectorType *y, const VectorType *x,
		INT SAFStart,
		const Configuration<MOType> & same) const;


	static MatrixType	getDiag(const Configuration<MOType> & same);

//---------------------------------------------------------------------	

static TwoIndexIntegralContainer *int2Container;
static FourIndexIntegralContainer *int4Container;
static const SymmetryContainer *iijj;		// pointer to iijj integrals

private:


static double	core;
static const MRFockMatrix<VectorType>	*mrFockMatrix;
static MatrixType	E0;
};





#endif
