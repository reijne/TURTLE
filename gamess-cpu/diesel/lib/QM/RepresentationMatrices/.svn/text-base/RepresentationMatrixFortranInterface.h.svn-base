//***********************************************************************
//
//	Name:			RepresentationMatrixFortranInterface.h
//
//	Description:	interface to the implementation in Fortran 
//					by Volker Pleﬂ
//					does the neccesary initializations and
//					calculates the representation matrices of the 
//					symmetric group approach (SGA) for one case
//					(matrices: "umat", "umat1", "emat", "hp5dar")
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.08.1996
//
//
//
//
//
//***********************************************************************

#ifndef __REPRESENTATIONMATRIXFORTRANINTERFACE_H
#define __REPRESENTATIONMATRIXFORTRANINTERFACE_H

#include "../../../config.h"

#include "../Math/SpinEigenFunctionDegeneration.h"
#include "RepresentationMatrix.h"



template <class RepMatType>
class RepresentationMatrixFortranInterface {
public:
	RepresentationMatrixFortranInterface(INT multiplicity);
	~RepresentationMatrixFortranInterface();

//----------------------------------------------------------------------	

/*	RepMatType *	getPRepMatCoulomb() const;
	RepMatType *	getPRepMatExchange() const;
	RepMatType *	getPRepMatP3() const;
	RepMatType *	getPRepMatP5() const;
*/
	
	INT	getSpinEigs(INT) const;
	INT	getMultiplicity() const;

//----------------------------------------------------------------------	

	void	calc(
		RepMatType *p1, RepMatType *p2,
		INT dK, INT P, INT R, INT qR, INT qL,
		INT kml, INT kmj,
		INT numberOfOpenShells
		);

//----------------------------------------------------------------------	

static SpinEigenFunctionDegeneration	*spinEigs;

private:
INT	matDim;
INT	maxOp;
double	*pUMat;
double	*pUMat1;
double	*pEMat;
//double	*php5dar;
};


/*
template <class RepMatType>
inline
RepMatType *	RepresentationMatrixFortranInterface<RepMatType>::getPRepMatCoulomb() const
{	return	pUMat;	}

template <class RepMatType>
inline
RepMatType *	RepresentationMatrixFortranInterface<RepMatType>::getPRepMatExchange() const
{	return	pUMat1;	}

template <class RepMatType>
inline
RepMatType *	RepresentationMatrixFortranInterface<RepMatType>::getPRepMatP3() const
{	return	pEMat;	}

//template <class RepMatType>
//inline
//RepMatType *	RepresentationMatrixFortranInterface<RepMatType>::getPRepMatP5() const
//{	return	php5dar;	}
*/

template <class RepMatType>
inline
INT	RepresentationMatrixFortranInterface<RepMatType>::getSpinEigs(INT i) const
{	return	(*spinEigs)(i);	}

template <class RepMatType>
inline
INT	RepresentationMatrixFortranInterface<RepMatType>::getMultiplicity() const
{	return	spinEigs->getMultiplicity();	}



#endif
