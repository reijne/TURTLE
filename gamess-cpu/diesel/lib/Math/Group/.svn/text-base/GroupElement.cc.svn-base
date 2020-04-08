//***********************************************************************
//
//	Name:			GroupElement.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27. Jan 1998
//
//***********************************************************************

#include <iostream>
#include <stdio.h>
#include <string>

#include "GroupElement.h"

using namespace std;



//*********************************************************************
//
// Implementation Groupelement
//

GroupElement::GroupElement()
{
	element = 0;
	group = (Group *) NULL;
}

GroupElement::GroupElement(Group & _group)
{
	element = 0;
	group = &_group;
}


GroupElement::~GroupElement()
{
	element = 0;
	group = (Group *) 0;
}



GroupElement & GroupElement::operator = (const GroupElement & a)
{
	group = a.getGroup();
	element = a.getElement();
	return *this;
}


GroupElement &GroupElement::operator *= (const GroupElement & a)
{
	if ( group != a.getGroup() )
	{	cerr << "GroupElement::operator *=: The factors belong to" <<
			" different groups.\n";
		return *this;
	}

	element = group->getProdukt(element, a.getElement());
	return *this;
}


GroupElement operator * 
	(const GroupElement & a, const GroupElement & b)
{
	if ( a.getGroup() != b.getGroup() )
	{	cerr << "GroupElement::operator *: The factors belong to" <<
			" different groups.\n";
		return a;
	}

GroupElement	g(*a.getGroup());
	g.setElement(a.getGroup()->
		getProdukt(a.getElement(), b.getElement()));
	return g;
}


GroupElement GroupElement::invert()
{
GroupElement	g(*group);

	g.setElement(group->getInv(element));
	return g;
}


ostream& operator<<(ostream& s, GroupElement & a)
{
	s << a.getElement();	
	return s;
}

/*
istream& operator>>(istream& s, GroupElement & a)
{
}
*/


