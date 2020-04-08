#include "../../../config.h"
//***********************************************************************
//
//	Name:			SymmetryContainer.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __SYMMETRYCONTAINER_H
#define __SYMMETRYCONTAINER_H


#include "IndexTranslation.h"
#include "BaseContainer.h"
#include "ExtPosContainer.h"
#include "../../Math/etc/iover2.h"

class SymmetryContainer : 
	public BaseContainer<ExtPosContainer> {
public:
	SymmetryContainer(INT i);
	~SymmetryContainer();

//--------------------------------------------------------------------	

	const NTupelContainer *	
		getNTupelContainer(const TwoElectronIntegralTriadeIndex &) const;

	const ExtPosContainer *	
		getExtPosContainer(const TwoElectronIntegralTriadeIndex &) const;

	//	the folowing operator does not perform any check if then integrals
	//	are actually contained in tree
	IntegralType	operator[] (const TwoElectronIntegralTriadeIndex &) const;


	//	all get/set methods return if the requested integral 
	//	is currently contained in tree
	INT	set(const TwoElectronIntegralTriadeIndex &, IntegralType);

	INT	get(const TwoElectronIntegralTriadeIndex &, IntegralType &) const;

	INT	get(const TwoElectronIntegralCbExIndex &,
		IntegralType &, IntegralType &) const;

//--------------------------------------------------------------------	


static MRMOs	mrmos;

};



inline
IntegralType	SymmetryContainer::operator[] (
	const TwoElectronIntegralTriadeIndex & ind) const
{
//	cout << "[]1" <<endl;
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		]->operator[](ind);

}


inline
INT	SymmetryContainer::get(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType &a) const
{
//	cout << "[]1" <<endl;
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		]->get(ind, a);

}


inline
INT	SymmetryContainer::set(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType a)
{
//	cout << "[]1" <<endl;
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		]->set(ind, a);

}


inline
INT	SymmetryContainer::get(
	const TwoElectronIntegralCbExIndex & ind,
	IntegralType & cb, IntegralType & ex) const
{
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		]->get(ind, cb, ex);

}


inline
const ExtPosContainer *	SymmetryContainer::getExtPosContainer(
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		];
}

inline
const NTupelContainer *	SymmetryContainer::getNTupelContainer(
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	iM1over2(
			iM1over2(mrmos.getIrRepP1(ind.getI())) + 
			mrmos.getIrRepP1(ind.getJ())
			) +
		iM1over2(mrmos.getIrRepP1(ind.getK())) + 
		mrmos.getIrRep(ind.getL())
		]->getNTupelContainer(ind);
}

#endif
