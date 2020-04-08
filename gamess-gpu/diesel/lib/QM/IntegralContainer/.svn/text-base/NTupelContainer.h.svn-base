#include "../../../config.h"
//***********************************************************************
//
//	Name:			NTupelContainer.h
//
//	Description:	base class for classes that contain
//					monades-, duades-, triade-integrals
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

#ifndef __NTUPELCONTAINER_H
#define __NTUPELCONTAINER_H

#include "IntegralType.h"

#include "IndexTranslation.h"
#include "../IntegralIndex/TwoElectronIntegralTriadeIndex.h"
#include "../IntegralIndex/TwoElectronIntegralCbExIndex.h"

class SharedMemory;

class NTupelContainer {
public:
	NTupelContainer(
		MRMOs &mrmos,
		INT i1,
		const IrRep *irreps,
		INT i2);
		
	~NTupelContainer();


//-----------------------------------------------------------------------

	//	the folowing operator does not perform any check if then integrals
	//	are actually contained in tree
	IntegralType operator[] (const TwoElectronIntegralTriadeIndex &) const;

	//	all get/set methods return if the requested integral 
	//	is currently contained in tree
	INT	set(const TwoElectronIntegralTriadeIndex &, IntegralType);

	INT	get(const TwoElectronIntegralTriadeIndex &, IntegralType &) const;

	INT	get(const TwoElectronIntegralCbExIndex &,
		IntegralType &, IntegralType &) const;

	
//-------------------------------------------------------------------------

	INT	isEmpty() const;

	void	check();

	INT	getDataSize() const;
	INT	allocate();

	void	setSharedMem(SharedMemory *sharedMem);

	INT	getNumberOfLeaves() const;	
	INT	getNumberOfIntegrals() const;	

//-----------------------------------------------------------------------



private:
	void	fillIdx(INT, INT *, INT);
	void	fillIdx(INT, INT *, INT *);

IntegralType	*p;				//	pointer to Integrals

INT				n;				//	number of contained integrals


INT				NTupel;			//	1=Monade, 2=Duade, 3=Triade

								//	for address calculation
const INT		*iIdx;			
const INT		*jIdx;			
const INT		*kIdx;			

INT				sub[4];			//	subtractors for address calculation

INT				map[3];			//	maps triade indices

const IndexTranslation *index;
SharedMemory	*sharedMem;		//	pointer to shared memory object
};


inline
INT	NTupelContainer::getDataSize() const
{	return n*3*sizeof(IntegralType);	}

/*
inline
void	NTupelContainer::get(
	TwoElectronIntegralTriadeIndex &ind, IntegralType & v)
{
INT	h = ind[indices];

	for ( INT i=0 ; i<indices ; i++ )
		h += mult[i] * ind[i];
		
	v = p[(h-extSub)*nTupel + ind.getM()];
}
*/


extern double hh;

inline
INT	NTupelContainer::isEmpty() const
{	return	!p;	}

inline
IntegralType	NTupelContainer::operator[] (
	const TwoElectronIntegralTriadeIndex & ind) const
{

//	printf("willi\n");
/*	printf("%x %d\n", p, n*NTupel);
	cout << "[]3" << p[0] <<endl;
	cout << "ind.getM()=" << ind.getM() <<", map[.])=" << map[ind.getM()] <<endl;
	cout << (iIdx[ind.getI()-sub[0]] +
				jIdx[ind.getJ()-sub[1]] +
				kIdx[ind.getK()-sub[2]] +
				ind.getL()-sub[3])*3 + map[ind.getM()] << endl;

	cout << iIdx[ind.getI()-sub[0]] << " " << 
				jIdx[ind.getJ()-sub[1]] << " " <<
				kIdx[ind.getK()-sub[2]] << " " <<
				ind.getL()-sub[3] << endl;

	cout << ind.getI()-sub[0] << " " << 
				ind.getJ()-sub[1] << " " <<
				ind.getK()-sub[2] << " " <<
				ind.getL()-sub[3] << endl;
*/

/*	cout << iIdx[ind[0]-sub[0]] << " " << 
				jIdx[ind[1]-sub[1]] << " " <<
				kIdx[ind[2]-sub[2]] << " " <<
				ind[3]-sub[3]  << endl;
	cout << ind[0]-sub[0] << " " << 
				ind[1]-sub[1] << " " <<
				ind[2]-sub[2] << " " <<
				ind[3]-sub[3]  << endl;*/
//	return hh;
	return p[	(iIdx[ind.getI()-sub[0]] +
				jIdx[ind.getJ()-sub[1]] +
				kIdx[ind.getK()-sub[2]] +
//				ind.getL()-sub[3])*NTupel + ind.getM()
				ind.getL()-sub[3])*3 + map[ind.getM()]
		];
}


inline
INT	NTupelContainer::set(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType a)
{
	if ( !p )
		return 0;
	p[	(iIdx[ind.getI()-sub[0]] +
				jIdx[ind.getJ()-sub[1]] +
				kIdx[ind.getK()-sub[2]] +
				ind.getL()-sub[3])*3 + map[ind.getM()]
		] = a;
		
	return 1;
}

inline
INT	NTupelContainer::get(
	const TwoElectronIntegralTriadeIndex & ind, IntegralType &a) const
{
	if ( !p )
		return 0;
	a = p[	(iIdx[ind.getI()-sub[0]] +
				jIdx[ind.getJ()-sub[1]] +
				kIdx[ind.getK()-sub[2]] +
				ind.getL()-sub[3])*3 + map[ind.getM()]
		];
	return 1;
}



inline
INT	NTupelContainer::get(
	const TwoElectronIntegralCbExIndex & ind,
	IntegralType &Cb, IntegralType &Ex) const
{
	if ( !p )
		return 0;
		
INT	h =	(iIdx[ind.getI()-sub[0]] +
			jIdx[ind.getJ()-sub[1]] +
			kIdx[ind.getK()-sub[2]] +
			ind.getL()-sub[3])*3;
			
	Cb = p[h + map[ind.getCoulombTriade()]];
	Ex = p[h + map[ind.getExchangeTriade()]];
	return 1;
}



inline
INT	NTupelContainer::getNumberOfIntegrals() const
{	return NTupel*n;	}


inline
INT	NTupelContainer::getNumberOfLeaves() const
{	return 1;	}




#endif
