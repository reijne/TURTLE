//***********************************************************************
//
//	Name:			Container.cc
//
//	Description:	abstract container
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.03.1997
//
//
//
//
//
//***********************************************************************




#include "Container.h"



template <class ContainedObjectType>
ostream& operator<<(ostream & s, const Container<ContainedObjectType> &c)
{
	s << c.getNumberOfElements() << std::endl;
	for ( ContainerIterator iter = c.first() ; !c.isLast(iter) ;
		c.next(iter) )
		if ( c.operator[](iter) )
			s << *c.operator[](iter);
	return s;
}

template <class ContainedObjectType>
void	Container<ContainedObjectType>::writeToStream(ostream & s) const
{
	s << getNumberOfNotNILElements() << std::endl;
	for ( ContainerIterator iter = first() ; !isLast(iter) ;
		next(iter) )
		if ( operator[](iter) )
			operator[](iter)->writeToStream(s);
}

template <class ContainedObjectType>
INT	Container<ContainedObjectType>::getNumberOfNotNILElements() const
{
INT	n = 0;
	for ( ContainerIterator iter = first() ; !isLast(iter) ;
		next(iter) )
		if ( operator[](iter) )
			n++;
	return n;
}

/*
template <class ContainedObjectType>
istream& operator>>(istream & s, Container<ContainedObjectType> &c)
{
INT	n;
	s >> n;
	for ( ContainerIterator iter = c.first() ; !c.isLast(iter) ;
		c.next(iter) )
	{
//		c.operator[](iter) = new ContainedObjectType&(n);
		s >> *c.operator[](iter);
	}
}
*/

#include "../Set/NExternalsSet.h"
#include "../Sel/NExternalsSel.h"
#include "../Diag/NExternalsDiag.h"
#include "../Set/InternalConfsSet.h"
#include "../Sel/InternalConfsSel.h"
#include "../Diag/InternalConfsDiag.h"
#include "../Set/TupelStructureSet.h"
#include "../Sel/TupelStructureSel.h"
#include "../Diag/TupelStructureDiag.h"
#include "../Set/extMOsSet.h"
#include "../Sel/extMOsSel.h"
#include "../Diag/extMOsDiag.h"
#include "../Base/extEntry.h"

template class Container<NExternalsSet>;
template class Container<NExternalsSel>;
template class Container<NExternalsDiag>;
template class Container<InternalConfsSet>;
template class Container<InternalConfsSel>;
template class Container<InternalConfsDiag>;
template class Container<TupelStructureSet>;
template class Container<TupelStructureSel>;
template class Container<TupelStructureDiag>;
template class Container<extMOsSet>;
template class Container<extMOsSel>;
template class Container<extMOsDiag>;
template class Container<extEntry>;


template ostream& operator << (ostream& s, const
    Container<NExternalsSet> & v);

template ostream& operator << (ostream& s, const
    Container<NExternalsSel> & v);

template ostream& operator << (ostream& s, const
    Container<NExternalsDiag> & v);


template ostream& operator << (ostream& s, const
    Container<InternalConfsSet> & v);

template ostream& operator << (ostream& s, const
    Container<InternalConfsSel> & v);

template ostream& operator << (ostream& s, const
    Container<InternalConfsDiag> & v);


template ostream& operator << (ostream& s, const
    Container<TupelStructureSet> & v);

template ostream& operator << (ostream& s, const
    Container<TupelStructureSel> & v);

template ostream& operator << (ostream& s, const
    Container<TupelStructureDiag> & v);


template ostream& operator << (ostream& s, const
    Container<extMOsSet> & v);

template ostream& operator << (ostream& s, const
    Container<extMOsSel> & v);

template ostream& operator << (ostream& s, const
    Container<extMOsDiag> & v);


template ostream& operator << (ostream& s, const
    Container<extEntry> & v);


//-------------------------------------------------------------------

/*
template class istream& operator >> (istream& s, 
    Container<InternalConfsSel> & v);

template class istream& operator >> (istream& s, 
    Container<InternalConfsDiag> & v);


template class istream& operator >> (istream& s, 
    Container<TupelStructureSel> & v);

template class istream& operator >> (istream& s, 
    Container<TupelStructureDiag> & v);


template class istream& operator >> (istream& s, 
    Container<extMOsSel> & v);

template class istream& operator >> (istream& s, 
    Container<extMOsDiag> & v);

template class istream& operator >> (istream& s, 
    Container<extEntry> & v);
*/


/*
istream& operator>>(istream & s, Container<InternalConfsSel> &c)
{
INT	n;
	s >> n;
	c.setNumberOfContainedObjects(n);
	for ( ContainerIterator iter = c.first() ; !c.isLast(iter) ;
		c.next(iter) )
	{
NExternalsSel *p = &c;
		c.operator[](iter) = new InternalConfsSel(p); //(&c);
		s >> *c.operator[](iter);
	}
}
*/
