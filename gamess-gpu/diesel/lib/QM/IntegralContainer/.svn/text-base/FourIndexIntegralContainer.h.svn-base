#include "../../../config.h"
//***********************************************************************
//
//	Name:			FourIndexIntegralContainer.h
//
//	Description:	root of tree containing four index integrals
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


#ifndef __FOURINDEXINTEGRALCONTAINER_H
#define __FOURINDEXINTEGRALCONTAINER_H


#include "../IntegralIndex/TwoElectronIntegralIndex.h"
#include "../IntegralIndex/TwoElectronIntegralTriadeIndex.h"
#include "../IntegralIndex/TwoElectronIntegralCbExIndex.h"

#include "IndexTranslation.h"
#include "BaseContainer.h"
#include "SymmetryContainer.h"

#include "IntegralType.h"

#include "../IO/Fortran/Fort31File.h"



class FourIndexIntegralContainer : public BaseContainer<SymmetryContainer> {
public:
	FourIndexIntegralContainer(const MRMOs & mrmos);
virtual	~FourIndexIntegralContainer();
	
//--------------------------------------------------------------------	


enum TIndDiffs {ijkl, ijkk, ijjk, ijjj, iijk, iijj, iiij, iiii};

	const NTupelContainer *	
		getNTupelContainer(const TwoElectronIntegralTriadeIndex &) const;

	const ExtPosContainer *	
		getExtPosContainer(const TwoElectronIntegralTriadeIndex &) const;

	const SymmetryContainer *	
		getSymmetryContainer(const TwoElectronIntegralTriadeIndex &) const;


	const SymmetryContainer *	operator[] (TIndDiffs i1) const;


	//	the folowing operators do not perform any check if then integrals
	//	are actually contained in tree
virtual	IntegralType	operator[] (const TwoElectronIntegralIndex<MOType> &) const;
virtual	IntegralType	operator[] (const TwoElectronIntegralTriadeIndex &) const;

	//	all get/set methods return if the requested integral 
	//	is currently contained in tree
	INT	set(const TwoElectronIntegralIndex<MOType> &, IntegralType);
	INT	set(const TwoElectronIntegralTriadeIndex &, IntegralType);

virtual	INT	get(const TwoElectronIntegralIndex<MOType> &, IntegralType &) const;
virtual	INT	get(const TwoElectronIntegralTriadeIndex &, IntegralType &) const;

virtual	INT	get(const TwoElectronIntegralCbExIndex &,
		IntegralType &, IntegralType &) const;



/*
	//	all get/set methods return if the requested integral 
	//	is currently contained in tree
	INT	get(
		TwoElectronIntegralIndex<MOType> &ind,
		IntegralType	&integral
		);
		
	INT	set(
		TwoElectronIntegralIndex<MOType> &ind,
		IntegralType	integral
		);


	INT	get(
		TwoElectronIntegralTriadeIndex &ind,
		IntegralType	&integral
		);

	INT	set(
		TwoElectronIntegralTriadeIndex &ind,
		IntegralType	integral
		);


	INT	get(
		TwoElectronIntegralCbExIndex &ind,
		IntegralType	&CbIntegral,
		IntegralType	&ExIntegral
		);

	INT	set(
		TwoElectronIntegralCbExIndex &ind,
		IntegralType	CbIntegral,
		IntegralType	ExIntegral
		);
*/	
	
//--------------------------------------------------------------------	

virtual	void	loadIntegrals(Fort31File f31);

//--------------------------------------------------------------------	
// 
// virtual	INT	getNumberOfTotalIntegrals() const;
// 
// virtual	INT	getNumberOfContainedIntegrals() const;


//--------------------------------------------------------------------	

	const MRMOs	*getMRMOs() const;

protected:
IndexTranslation	*index;
MRMOs	mrmos;
};




inline
const MRMOs	*FourIndexIntegralContainer::getMRMOs() const
{
	return	&mrmos;
}

inline
IntegralType	FourIndexIntegralContainer::operator[] (
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->operator[](ind);
}

inline
IntegralType	FourIndexIntegralContainer::operator[] (
	const TwoElectronIntegralIndex<MOType> & ind) const
{
TwoElectronIntegralTriadeIndex	triadeInd(ind);

	return p[	triadeInd.getK()==triadeInd.getL() |
		((triadeInd.getJ()==triadeInd.getK()) << 1) |
		((triadeInd.getI()==triadeInd.getJ()) << 2)
		]->operator[](triadeInd);
}

inline
INT	FourIndexIntegralContainer::get(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType &a) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->get(ind, a);
}

inline
INT	FourIndexIntegralContainer::get(
	const TwoElectronIntegralIndex<MOType> & ind, IntegralType &a) const
{
TwoElectronIntegralTriadeIndex	triadeInd(ind);

	return p[	triadeInd.getK()==triadeInd.getL() |
		((triadeInd.getJ()==triadeInd.getK()) << 1) |
		((triadeInd.getI()==triadeInd.getJ()) << 2)
		]->get(triadeInd, a);
}

inline
INT	FourIndexIntegralContainer::set(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType a)
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->set(ind, a);
}

inline
INT	FourIndexIntegralContainer::set(
	const TwoElectronIntegralIndex<MOType> & ind, IntegralType a)
{
TwoElectronIntegralTriadeIndex	triadeInd(ind);

	return p[	triadeInd.getK()==triadeInd.getL() |
		((triadeInd.getJ()==triadeInd.getK()) << 1) |
		((triadeInd.getI()==triadeInd.getJ()) << 2)
		]->set(triadeInd, a);
}

inline
const SymmetryContainer *	FourIndexIntegralContainer::operator[]
	(TIndDiffs i1) const
{	return p[i1];	}


inline
INT	FourIndexIntegralContainer::get(
	const TwoElectronIntegralCbExIndex & ind,
	IntegralType & cb, IntegralType & ex) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->get(ind, cb ,ex);
}


inline
const SymmetryContainer *	FourIndexIntegralContainer::getSymmetryContainer(
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		];
}

inline
const ExtPosContainer *	FourIndexIntegralContainer::getExtPosContainer(
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->getExtPosContainer(ind);
}

inline
const NTupelContainer *	FourIndexIntegralContainer::getNTupelContainer(
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	ind.getK()==ind.getL() |
		((ind.getJ()==ind.getK()) << 1) |
		((ind.getI()==ind.getJ()) << 2)
		]->getNTupelContainer(ind);
}


#endif
