#include "../../../config.h"
//***********************************************************************
//
//	Name:			ExtPosContainer.h
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

#ifndef __EXTPOSCONTAINER_H
#define __EXTPOSCONTAINER_H

#include "IndexTranslation.h"
#include "BaseContainer.h"
#include "NTupelContainer.h"

#include "../MO/MRMOs.h"

class ExtPosContainer : public BaseContainer<NTupelContainer> {
public:
	ExtPosContainer(INT i1, const IrRep *irreps);

	~ExtPosContainer();

//--------------------------------------------------------------------	

	const NTupelContainer *	
		getNTupelContainer(const TwoElectronIntegralTriadeIndex &) const;

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

private:
};

inline
IntegralType 	ExtPosContainer::operator[] (
	const TwoElectronIntegralTriadeIndex & ind) const
{
/*	cout << "[]2" <<endl;	
	
	cout << (mrmos.isInternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)) << endl;
		
	printf("%x\n", p[mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)]);
*/		
		
	return p[	mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)
		]->operator[](ind);

}


inline
INT	ExtPosContainer::get(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType &a) const
{
	return p[	mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)
		]->get(ind, a);
}


inline
INT	ExtPosContainer::set(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType a)
{
	return p[	mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)
		]->set(ind, a);
}


inline
INT	ExtPosContainer::get(
	const TwoElectronIntegralCbExIndex & ind,
	IntegralType & cb, IntegralType & ex) const
{
	return p[	mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)
		]->get(ind, cb, ex);

}


inline
const NTupelContainer *	ExtPosContainer::getNTupelContainer (
	const TwoElectronIntegralTriadeIndex & ind) const
{
	return p[	mrmos.isExternal(ind.getL()) |
		(mrmos.isExternal(ind.getK()) << 1) |
		(mrmos.isExternal(ind.getJ()) << 2) |
		(mrmos.isExternal(ind.getI()) << 3)
		];
}


#endif
