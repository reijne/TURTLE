//***********************************************************************
//
//	Name:			InternalConfsSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.05.1997
//
//
//
//
//
//***********************************************************************


#ifndef __InternalConfsSet_H
#define __InternalConfsSet_H

#include "../../../../config.h"


#include "../Base/MRTreeBase.h"
#include "../Base/InternalConfsBase.h"

//#include "TupelStructureSet.h"

#include "../Container/SetContainerTree.h"

template <class TMOType> class Configuration;
class NExternalsSet;

class InternalConfsSet : 
	public InternalConfsBase,
	public SetContainerTree<NExternalsSet, TupelStructureSet>,
	public MRTreeBase<NExternalsSet, TupelStructureSet> {
public:
	InternalConfsSet(INT nExt, NExternalsSet *parent) :
                Tree<NExternalsSet>(parent),
		InternalConfsBase(nExt)  {   generated = 0;}
//	InternalConfsSet(NExternalsSet *parent);
	InternalConfsSet(istream &s);
	~InternalConfsSet();
	

	// "virtual" constructor:
//	virtual InternalConfsSet * new_Set();

	void	copy(InternalConfsSet *, INT deep = 0);
	virtual SetContainerTree<NExternalsSet, TupelStructureSet> *
		clone(INT deep = 0);
	
	
//	InternalConfsSet(const InternalConfsSet &);
	InternalConfsSet & operator = (const InternalConfsSet &);
	
	
	//----------------------------------------------------------------	

	void	add(const Configuration<MOType> &internal,
		const Configuration<MOType> &external);

	void	add(TupelStructureSet *);

	//----------------------------------------------------------------	

	INT	getNumberOfGenerated() const;

	INT checkSameRefs(const InternalConfsSet &) const;

	void	cutTreeLevel2();

	//--------------------------------------------------------------------------

	INT operator == (const InternalConfsSet &) const;
	INT operator != (const InternalConfsSet &) const;
	INT operator <= (const InternalConfsSet &) const;

	//----------------------------------------------------------------
	
//	void	operator |= (InternalConfsSet& b);

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const InternalConfsSet &);
	friend istream& operator>>(istream & s, InternalConfsSet &);

	//----------------------------------------------------------------	

private:
INT	generated;		//	number of generated configurations
};



inline
INT	InternalConfsSet::getNumberOfGenerated() const
{	return	generated;	}

inline
void	InternalConfsSet::add(TupelStructureSet *t)
{
	SetContainerTree<NExternalsSet, TupelStructureSet>::add(t);
}

/*
inline
InternalConfsSet::InternalConfsSet(const InternalConfsSet &a) :
	SetContainer<TupelStructureSet>(),	
	InternalConfsBase(a)
{			cout << "!!!!:" << *(operator [] (first())) << endl;
	generated = a.getNumberOfGenerated();	}


*/
inline
InternalConfsSet & InternalConfsSet::operator = (const InternalConfsSet &a)
{
//			cout << "!!!!InternalConfsSet::operator =:" << *(operator [] (first())) << endl;
	generated = a.getNumberOfGenerated();
	(InternalConfsBase &) *this	= a;
	return *this;
}

/*
inline 
InternalConfsSet *	InternalConfsSet::
	new_Set()
{	return new InternalConfsSet(nExt, NULL);	}
*/

inline 
void	InternalConfsSet::copy(InternalConfsSet *c, INT deep)
{
//	cout << "InternalConfsSet::copy:" << *c << endl;
	*this = *c;
//	cout << *this << endl;
	if ( deep )
		SetContainerTree<NExternalsSet, TupelStructureSet>::clone(1);
}


inline 
SetContainerTree<NExternalsSet, TupelStructureSet> *	InternalConfsSet::
	clone(INT deep)
{	
	InternalConfsSet * c = new InternalConfsSet(0, NULL);
//	cout << "InternalConfsSet::clone(INT deep)" << endl;
	c->copy(this, deep);
//	cout << "after copy" << endl;
	return c;
}


inline
INT InternalConfsSet::operator == (const InternalConfsSet &a) const
{	return InternalConfsBase::operator == (a);	}

inline
INT InternalConfsSet::operator != (const InternalConfsSet &a) const
{	return InternalConfsBase::operator != (a);	}

inline
INT InternalConfsSet::operator <= (const InternalConfsSet &a) const
{	return InternalConfsBase::operator <= (a);	}

#endif
