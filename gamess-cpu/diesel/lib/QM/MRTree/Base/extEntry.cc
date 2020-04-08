//***********************************************************************
//
//	Name:			extEntry.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.03.1997
//
//
//
//
//
//***********************************************************************

#include "extEntry.h"

#include "../../MO/MOMapping.h"

#include <iomanip>

#include <stdlib.h>


#include "extMOsBase.h"
#include "TupelStructureBase.h"

using namespace std;

//	Bug in GNU C++ 2.7.2.3 on AIX: Classes with virtual functions must have
//	at least one virtual function defined non-inlined
extEntry::~extEntry()
{
	if ( mo )
		delete mo;	
}





//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//	Attention:
//
//	MO sorting scheme due to MO mapping will not function properly
//	for more than 2 external MOs because it will mix open and closed
//	shells as there is no information about # open and closed shells
//	in this object.
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



static int	cmp(const void *_p1, const void *_p2)
{
const MOType *p1 = (const MOType *) _p1;
const MOType *p2 = (const MOType *) _p2;
	if ( *p1<*p2 )
		return -1;
	if ( *p1>*p2 )
		return 1;
	return 0;
}
		

extEntry::extEntry(istream &s) :
	RootEnergies(s), Tree<extMOsBase>(NULL)
{
	s >> n;
	mo = NULL;
	if ( !n )
		return;
		
	mo = new MOType[n];
	for ( INT i=0 ; i<n ; i++ )
	{
		s >> mo[i];
		mo[i] = moMapping.getContinuous(mo[i]);
	}
	qsort(mo, n, sizeof(MOType), cmp);

	s >> SelectionFlags;

//	cout << "::::::::::::::::" << endl;
//	cout << *this;
//	cout << ";;;;;;;;;;;;;;;;" << endl;
	
}


void	extEntry::writeToStream(ostream & s) const
{
	RootEnergies::writeToStream(s);
	s << n << endl;

MOType* m = new MOType[n];
	for ( INT i=0 ; i<n ; i++ )
		m[i] = moMapping.getReal(mo[i]);

	qsort(m, n, sizeof(MOType), cmp);
		
	for ( INT i=0 ; i<n ; i++ )
		s << m[i] << " ";

	if ( n )
	{
		s << endl;
		s << SelectionFlags << endl;
	}
	delete[] m;
}


ostream& operator<<(ostream & s, const extEntry &e)
{
	s << e.n << endl;

MOType* m = new MOType[e.n];
	for ( INT i=0 ; i<e.n ; i++ )
		m[i] = moMapping.getReal(e.mo[i]);

	qsort(m, e.n, sizeof(MOType), cmp);
		
	for ( INT i=0 ; i<e.n ; i++ )
		s << m[i] << " ";
	s << endl;
	
	s << e.nRoots << endl;
	for ( INT i=0 ; i<e.nRoots ; i++ )
		s << setiosflags(ios::scientific) << e.e[i] << " ";
	s << endl;

	s << e.SelectionFlags << endl;
	delete[] m;
	return s;
}


istream& operator>>(istream & s, extEntry &e)
{
	s >> e.n;
	for ( INT i=0 ; i<e.n ; i++ )
	{
		s >> e.mo[i];
		e.mo[i] = moMapping.getContinuous(e.mo[i]);
	}
	qsort(e.mo, e.n, sizeof(MOType), cmp);

	s >> e.nRoots;
	for ( INT i=0 ; i<e.nRoots ; i++ )
		s >> setiosflags(ios::scientific) >> e.e[i];

	s >> e.SelectionFlags;
	return s;
}


extEntry &	extEntry::operator |= (extEntry & b)
{
	return *this;
}

extEntry &	extEntry::operator &= (extEntry & b)
{
	return *this;
}

extEntry &	extEntry::operator -= (extEntry & b)
{
	return *this;
}



Configuration<MOType>	extEntry::getConfiguration() const
{
Configuration<MOType>	conf = *(getParent()->getParent());

//	cout << "n= " << n << ", this = " << this << endl;
INT	j = 0;
	for ( INT i=0 ; i<getParent()->getNumberOfOpenMOs() ; i++ )
		conf.create(mo[j++]);

	for ( INT i=0 ; i<getParent()->getNumberOfClosedMOs() ; i++ )
	{
		conf.create(mo[j]);
		conf.create(mo[j++]);
	}

	return conf;
}




#include "../../../Container/SLList.cc"
template class SLList<extEntry *>;


#define AVLSET_CLONE
#include "../../../Container/AVLSet.cc"
template class AVLSet<extEntry *>;

#include "../../../Container/Set.cc"
template class Set<extEntry *>;
