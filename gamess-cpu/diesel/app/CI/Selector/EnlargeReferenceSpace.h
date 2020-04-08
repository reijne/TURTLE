//***********************************************************************
//
//	Name:			EnlargeReferenceSpace.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			01.04.1998
//
//
//
//
//
//***********************************************************************



#ifndef __EnlargeReferenceSpace_h
#define __EnlargeReferenceSpace_h

#include "../../../config.h"

#include "../../../lib/QM/MRTree/EnergyType.h"
#include "../../../lib/QM/RepresentationMatrices/DataTypes.h"


#include "../../../lib/QM/MO/MOType.h"
#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/RepresentationMatrices/HMatElements.h"

class MRConfInput;
class NExternalsSet;
/* FD class ostream; */
class BinomialCoefficient;
template <class MOType> class	TableCase;
class TableKey;
template <class TableKey, class RepresentationMatrices> class	VarSizeReadOnlyCache;
class SpinEigenFunctionDegeneration;



template <class RepMatType, class VectorType>
class EnlargeReferenceSpace {
public:
	EnlargeReferenceSpace(const MRConfInput &mrinp);
	~EnlargeReferenceSpace();

	
	const EnergyType *	calcEnlargementEnergy(const Configuration<MOType> &);

	RootType	getRootEnergy(INT i) const;


        template <class RT,class VT>
	friend ostream & operator << (ostream &s, const EnlargeReferenceSpace<RT, VT> &e);

private:
	EnlargeReferenceSpace(const EnlargeReferenceSpace<RepMatType, VectorType> &);
	EnlargeReferenceSpace & operator = (const EnlargeReferenceSpace<RepMatType, VectorType> &);


NExternalsSet	*PT0Wave;
//ConfigurationSet confSet;			//	set of reference configurations
INT	dim;							//	reference space dimension
RootType	*refEnergies;			//	energies from reference space
INT	nRoots;							//	number of roots
INT	*roots;							//	field of nth root
EnergyType	*enlargedEnergies;		//	energies after appending additional configuration
RepMatType	*pRefMat;				//	reference matrix
SpinEigenFunctionDegeneration	*spinEigs; 
BinomialCoefficient	*binom;
TableCase<MOType>	*tablecase;
VarSizeReadOnlyCache<TableKey, HMatElements<RepMatType, VectorType> >	*cache;
};




template <class RepMatType, class VectorType>
inline
RootType	EnlargeReferenceSpace<RepMatType, VectorType>::getRootEnergy(INT i) const
{	return refEnergies[i];	}


#endif
