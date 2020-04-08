#include "../../../config.h"
//***********************************************************************
//
//	Name:			TwoIndexIntegralContainer.h
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

#ifndef __TWOINDEXINTEGRALCONTAINER_H
#define __TWOINDEXINTEGRALCONTAINER_H


#include "../IntegralIndex/OneElectronIntegralIndex.h"

#include "IntegralType.h"
#include "../MO/MRMOs.h"

#include "../IO/Fortran/Fort31File.h"
#include <iostream>

class TwoIndexIntegralContainer {
public:
	TwoIndexIntegralContainer(const MRMOs & mrmos);
	virtual ~TwoIndexIntegralContainer();
	
//--------------------------------------------------------------------	



	OneElectronIntegralType &	operator[] (const OneElectronIntegralIndex &);


	
//--------------------------------------------------------------------	

	virtual double	loadIntegrals(Fort31File f31);

//--------------------------------------------------------------------	

	INT	getNumberOfTotalIntegrals() const;

	INT	getNumberOfContainedIntegrals() const;

	double getCore() const
	{
		//cout << "*** core-Energie= " << core << endl;
		return core;
	} 

//--------------------------------------------------------------------	


protected:
OneElectronIntegralType	*p;				//	pointer to Integrals
INT						n;				//	number of contained integrals
INT						*tab;			//	contains sum_{i=1}^n i
const MRMOs	*mrmos;
double	core;
};



inline
OneElectronIntegralType &	TwoIndexIntegralContainer::operator[] (
	const OneElectronIntegralIndex & ind)
{
	if ( ind.getI()>=ind.getJ() )
	{
		return p[tab[ind.getI()-1] + ind.getJ() - 1];
	}

	else
	{
		return p[tab[ind.getJ()-1] + ind.getI() - 1];
	}
}



#endif
