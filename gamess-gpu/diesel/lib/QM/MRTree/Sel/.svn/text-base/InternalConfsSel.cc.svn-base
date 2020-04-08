//***********************************************************************
//
//	Name:			InternalConfsSel.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************

#include "NExternalsSel.h"
#include "InternalConfsSel.h"
#include "TupelStructureSel.h"

#include <stdio.h>


/*
InternalConfsSel::InternalConfsSel(NExternalsSel *parent) :
	Tree<NExternalsSel>(parent)
{
	generated = 0;
}
*/

InternalConfsSel::InternalConfsSel(istream &s) :
        Tree<NExternalsSel>(NULL),
	SetContainer<TupelStructureSel>(s),
        InternalConfsBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
	generated = 0;
}


InternalConfsSel::~InternalConfsSel()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		delete operator [] (iter);
}


void	InternalConfsSel::cutTreeLevel2()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
	{
		if ( operator [] (iter) )
		{
			if ( operator [] (iter)->getNumberOfLeaves()==0 )
			{
				delete operator [] (iter);
				operator [] (iter) = NULL;
			}
		}
	}
}


void	InternalConfsSel::addExtEntry(
	const Configuration<MOType> & conf, const RootEnergies & r)
{
Configuration<MOType>	internal;
Configuration<MOType>	external;

//	cout << "HALLO: InternalConfsSel::addExtEntry" << endl;
	conf.split(getParent()->getMRMOs(), internal, external);

	for ( ContainerIterator iter=first() ; !isLast(iter) ; next(iter) )
		if ( internal == *(*this)[iter] )
			(*this)[iter]->addExtEntry(external, r);
}

const extEntry *	InternalConfsSel::getExtEntry(
	Configuration<MOType> conf) const
{
Configuration<MOType>	internal;
Configuration<MOType>	external;

//	cout << "HALLO: InternalConfsSel::getExtEntry" << endl;
	conf.split(getParent()->getMRMOs(), internal, external);

	for ( ContainerIterator iter=first() ; !isLast(iter) ; next(iter) )
		if ( internal == *(*this)[iter] )
			return (*this)[iter]->getExtEntry(external);


	return NULL;
}


void	InternalConfsSel::threshCut(EnergyType Threshold, INT divideByCSFs)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			if ( !operator [] (iter)->isReference() )
				operator [] (iter)->threshCut(Threshold, divideByCSFs);
}


/*
void	InternalConfsSel::setNumberOfConfigurations(INT n)
{
	number = n;
	mainConfiguration = new TupelStructureSel[number];
}
*/

void	InternalConfsSel::writeToStream(ostream & s) const
{
	Container<TupelStructureSel>::writeToStream(s);
	InternalConfsBase::writeToStream(s);
}

ostream& operator<<(ostream & s, const InternalConfsSel &base)
{
	s << ((InternalConfsBase &) base);
	s << ((Container<TupelStructureSel> &) base);
	return s;	
}

/*
istream& operator>>(istream & s, InternalConfsSel &base)
{
	s >> ((InternalConfsBase &) base);
	s >> ((Container<TupelStructureSel> &) base);
	return s;	
}
*/
