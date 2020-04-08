//***********************************************************************
//
//	Name:			PTReferenceSpace.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.04.1997
//
//
//
//
//
//***********************************************************************



#ifndef __PTReferenceSpace_h
#define __PTReferenceSpace_h

#include "../../../config.h"

#include "../../../app/CI/Selector/MRConfInput.h"
#include "../../Math/etc/Histogram.h"
#include "../MRTree/EnergyType.h"

#include "DataTypes.h"


template <class DensMatType> class MRFockMatrix;

class	NExternalsSet;
class	PTSum;

template <class MatrixType, class VectorType> class CICalculation;

template <class MatrixType, class VectorType>
class PTReferenceSpace {
public:
	PTReferenceSpace(MRConfInput mrinp,
		const CICalculation<MatrixType, VectorType> &);
	~PTReferenceSpace();
	
	
	INT	getNumberOfRoots() const;
	RootType	getRoot(INT rootNr) const;
	const CoefType *getCoefP(INT confNr, INT rootNr) const;
	CoefType getCoef(INT confNr, INT rootNr, INT SAFNr) const;

	const NExternalsSet	*getPT0Wave() const;
	
	
	Histogram<EnergyType>	*getHistogram(INT i) const;
	Histogram<EnergyType>	**getHistogramP() const;
	PTSum *getPTTotalSum() const;
	PTSum *getPTSum(INT i) const;
	PTSum **getPTSumP() const;
	
	ConfigurationSet getRefSet() const;
	
private:
INT	nRoots;					// number of roots
RootType	*roots;			// pointer to array of roots
CoefType	***coefs;		// pointer to array of coefs
Histogram<EnergyType>	**hist;
							// pointer to array of histograms

PTSum	*ptTotalSum;		// selection information (selection procedure
							// or'ed with respect to different roots)
							// concerning several thresholds

PTSum	**ptSum;			// pointer to array of selection information
							// concerning several thresholds

NExternalsSet	*PT0Wave;	// pointer to tree of configurations

VectorType	**densMats;		// pointer to reference density matrices

MRFockMatrix<VectorType>	*mrFockMatrix;


ConfigurationSet	refSet;
};


template <class MatrixType, class VectorType>
inline
INT	PTReferenceSpace<MatrixType, VectorType>::getNumberOfRoots() const
{	return nRoots;	}

template <class MatrixType, class VectorType>
inline
RootType	PTReferenceSpace<MatrixType, VectorType>::getRoot(INT i) const
{	return roots[i];	}

template <class MatrixType, class VectorType>
inline
const NExternalsSet *	PTReferenceSpace<MatrixType, VectorType>::getPT0Wave() const
{	return PT0Wave;	}

template <class MatrixType, class VectorType>
inline
Histogram<EnergyType> *	PTReferenceSpace<MatrixType, VectorType>::getHistogram(INT i) const
{	return hist[i];	}

template <class MatrixType, class VectorType>
inline
Histogram<EnergyType> **	PTReferenceSpace<MatrixType, VectorType>::getHistogramP() const
{	return hist;	}

template <class MatrixType, class VectorType>
inline
PTSum *	PTReferenceSpace<MatrixType, VectorType>::getPTSum(INT i) const
{	return ptSum[i];	}

template <class MatrixType, class VectorType>
inline
PTSum **	PTReferenceSpace<MatrixType, VectorType>::getPTSumP() const
{	return ptSum;	}

template <class MatrixType, class VectorType>
inline
PTSum *	PTReferenceSpace<MatrixType, VectorType>::getPTTotalSum() const
{	return ptTotalSum;	}


template <class MatrixType, class VectorType>
inline
const CoefType *PTReferenceSpace<MatrixType, VectorType>::getCoefP(INT confNr, INT rootNr) const
{	return coefs[confNr][rootNr];	}

template <class MatrixType, class VectorType>
inline
CoefType PTReferenceSpace<MatrixType, VectorType>::getCoef(INT confNr, INT rootNr, INT SAFNr) const
{	return coefs[confNr][rootNr][SAFNr];	}

template <class MatrixType, class VectorType>
inline
ConfigurationSet PTReferenceSpace<MatrixType, VectorType>::getRefSet() const
{	return refSet;	}

#endif

