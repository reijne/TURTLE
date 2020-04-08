#include "../../../config.h"
//***********************************************************************
//
//	Name:			MRFockMatrix.h
//
//	Description:	generalized fock matrix to multi reference case
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			15.09.1998
//
//	Literature:		-	Robert B. Murphy, Richard P. Messmer:
//						"Generalized Moller-Plesset perturbation theory
//						applied to general MCSCF reference wave functions",
//						Chem. Phys. Lett., 183/5, (1991), p. 443
//					-	Krzysztof Wolinski, Peter Pulay:
//						"Generalized Moller-Plesset perturbation theory:
//						Second order results for two-configuration, open-shell
//						excited singlet, and doublet wave functions",
//						J. Chem. Phys. 90 (7), 1989, 3647
//					-	Krzysztof Wolinski, Harrell L. Sellers, Peter Pulay:
//						"Consistent generalization of the Moller-Plesset
//						partitioning to open-shell and multiconfigurational SCF
//						reference states in many-body perturbation theory",
//						Chem. Phys. Lett. 140/3, 1987, 225
//
//
//
//
//***********************************************************************

#ifndef __MRFockMatrix_h
#define __MRFockMatrix_h


#include "IntegralType.h"
#include "../MO/MOType.h"


class FourIndexIntegralContainer;
class TwoIndexIntegralContainer;
class NExternalsDiag;

template <class MatrixType, class VectorType> class CICalculation;

template <class DensMatType>
class MRFockMatrix {
public:
	// initialize static calling tables

	//	for hamilton matrix multiplication
	template <class MatrixType> MRFockMatrix(
		const CICalculation<MatrixType, DensMatType> *,
		const NExternalsDiag *,
		const DensMatType *alpha,
		FourIndexIntegralContainer *,
		TwoIndexIntegralContainer *);

	~MRFockMatrix();


	IntegralType	operator () (MOType p, MOType q) const;
	

private:
	
	
INT	maxMO;

								// -----
                             	//  \       #root
IntegralType	*mrFockMatrix;	//   >     D      [(pq|rs) - (pr|qs)/2]
								//  /       rs
								// -----
								//  r,s
								//
								//  r, s: internal occupied orbitals
								//
								//  D   : reference density
								//   rs
								//
								//  p, q: running indices in mrFockMatrix
};



template <class DensMatType>
inline
IntegralType	MRFockMatrix<DensMatType>::operator ()
	(MOType p, MOType q) const
{
	return mrFockMatrix[(p-1)*maxMO + (q-1)];
}


#include "MRFockMatrix.cch"


#endif
