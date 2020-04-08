//***********************************************************************
//
//	Name:			extMOsBase.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			11.03.1997
//
//
//
//
//
//***********************************************************************

#include "NExternalsBase.h"
#include "InternalConfsBase.h"
#include "TupelStructureBase.h"
#include "extMOsBase.h"

#include <stdlib.h>
#include <string>

extMOsBase::extMOsBase(istream &s)
{
//	cout << "extMOsBase(istream &s)" << endl;
	s >> total;
	s >> open;
	s >> SAFStart;
	s >> SAFInc;
	s >> Syms;
	inSymN = new INT [Syms];
	for ( INT i=0 ; i<Syms ; i++ )
		s >> inSymN[i];
	s >> maxMO;
//	cout << "OOOOOOOOOOOOOOOOOOOOO" << endl;
//	cout << *this << endl;
//	cout << ",,,,,,,,,,,,,,,,,,,,," << endl;
}

extMOsBase::extMOsBase()
{
	inSymN = NULL;
	total = open = SAFStart = SAFInc = Syms = maxMO = 0;
}

extMOsBase::extMOsBase(
	INT _Syms,
	INT _open,
	INT _closed)
{
	Syms = _Syms;
	total = _open + _closed;
	open = _open;
	inSymN = new INT [Syms];
	if ( Syms==1 )
		inSymN[0] = 0;
        SAFStart = SAFInc = maxMO = 0;
}

extMOsBase::~extMOsBase()
{
	if ( inSymN )
		delete[] inSymN;
        inSymN = NULL;
}


extMOsBase::extMOsBase(const extMOsBase &ext)
{
	total = ext.getNumberOfOpenMOs() + ext.getNumberOfClosedMOs();
	open = ext.getNumberOfOpenMOs();
	SAFStart = ext.getSAFStart();
	SAFInc = ext.getSAFInc();
	maxMO = ext.getMaxMO();
	Syms = ext.getNumberOfSymmetries();
        if (ext.inSymN) 
        {
	   inSymN = new INT [Syms];
	   for ( INT i=0 ; i<Syms ; i++ )
		inSymN[i] = ext.getInSymN(i);
        }
        else
           inSymN = NULL;

}


extMOsBase & extMOsBase::operator = (const extMOsBase &ext)
{
	total = ext.getNumberOfOpenMOs() + ext.getNumberOfClosedMOs();
	open = ext.getNumberOfOpenMOs();
	SAFStart = ext.getSAFStart();
	SAFInc = ext.getSAFInc();
	maxMO = ext.getMaxMO();
	Syms = ext.getNumberOfSymmetries();
	if ( inSymN )
		delete inSymN;
        if (ext.inSymN) 
        {
	   inSymN = new INT [Syms];
	   for ( INT i=0 ; i<Syms ; i++ )
		 inSymN[i] = ext.getInSymN(i);
        }
        else
           inSymN = NULL;
	return *this;
}



const extEntry *	extMOsBase::getExtEntry(
	Configuration<MOType>  external) const
{
INT	n = external.getNumberOfOpenShells() + external.getNumberOfClosedShells();
MOType* mo = new MOType[n];
INT	j = 0;
//	cout << "HALLO: extMOsBase::getExtEntry" << endl;
//	cout << ((Configuration<MOType>) *getParent()) << " ::" << getNumberOfElements()<< endl;
	for ( INT i=0 ; i<external.getNumberOfOpenShells() ; i++ )
		mo[j++] = external.getOpenShell(i);
	for ( INT i=0 ; i<external.getNumberOfClosedShells() ; i++ )
		mo[j++] = external.getClosedShell(i);

	for ( ContainerIterator iter=first() ; !isLast(iter) ; next(iter) )
	{
//		cout << (*this)[iter]->getConfiguration() << endl;
		if ( !memcmp(mo, (*this)[iter]->getMOP(), n*sizeof(MOType)) )
			return (*this)[iter];
	}

	delete[] mo;
	return NULL;
}


void	extMOsBase::writeToStream(ostream & s) const
{
	s << total << endl;
	s << open << endl;
	s << SAFStart << endl;
	s << SAFInc << endl;
	s << Syms << endl;
	for ( INT i=0 ; i<Syms ; i++ )
		s << inSymN[i] << " ";
	s << endl;
	s << maxMO << endl;
}


ostream& operator<<(ostream & s, const extMOsBase & mo)
{
	s << mo.total << endl;
	s << mo.open << endl;
	s << mo.SAFStart << endl;
	s << mo.SAFInc << endl;
	s << mo.Syms << endl;
	for ( INT i=0 ; i<mo.Syms ; i++ )
		s << mo.inSymN[i] << " ";
	s << endl;
	s << mo.maxMO << endl;
/*	for ( INT i=0 ; i<mo.getNumberOfConfs() ; i++ )
	{	for ( INT j=0 ; j<mo.getNumberOfOpenMOs() ; j++ )
			s << mo.getOpenMO(i, j) << " ";
		if ( mo.getNumberOfOpenMOs() && mo.getNumberOfClosedMOs() )
			s << "| ";
		for ( INT j=0 ; j<mo.getNumberOfClosedMOs() ; j++ )
			s << mo.getClosedMO(i, j) << " ";
		s << endl;
	}			
*/
	return s;
}

istream& operator>>(istream & s, extMOsBase & mo)
{
	s >> mo.total;
	s >> mo.open;
	s >> mo.SAFStart;
	s >> mo.SAFInc;
	s >> mo.Syms;
	for ( INT i=0 ; i<mo.Syms ; i++ )
		s >> mo.inSymN[i];
	s >> mo.maxMO;
	return s;
}

