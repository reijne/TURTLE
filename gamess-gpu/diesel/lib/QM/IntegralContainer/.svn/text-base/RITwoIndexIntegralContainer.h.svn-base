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

#ifndef __RITWOINDEXINTEGRALCONTAINER_H
#define __RITWOINDEXINTEGRALCONTAINER_H


#include "TwoIndexIntegralContainer.h"
#include <string>
#include <iostream>

#include "../IntegralIndex/OneElectronIntegralIndex.h"

#include "IntegralType.h"
#include "../MO/MRMOs.h"

#include "../IO/Fortran/Fort31File.h"

class RITwoIndexIntegralContainer:public TwoIndexIntegralContainer {
public:
	RITwoIndexIntegralContainer(const MRMOs & mrmos);
	~RITwoIndexIntegralContainer();
	

	double	loadIntegrals(Fort31File f31);


//--------------------------------------------------------------------	
// geerbte Daten aus TwoIndexIntegralContainer
// private:
// OneElectronIntegralType	*p;				//	pointer to Integrals
// INT						n;				//	number of contained integrals
// INT						*tab;			//	contains sum_{i=1}^n i
// const MRMOs	*mrmos;
// double	core;
//--------------------------------------------------------------------	
private:
string 		oneintFilename;
};



#endif
