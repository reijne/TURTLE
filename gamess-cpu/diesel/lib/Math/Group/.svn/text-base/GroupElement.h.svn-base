//***********************************************************************
//
//	Name:			GroupElement.h
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




#ifndef _GROUPELEMENT_H
#define _GROUPELEMENT_H

#include "../../../config.h"

#include "Group.h"



class GroupElement {
public:
	GroupElement();
	GroupElement(Group &g);
	~GroupElement();
	
	
	GroupElement & operator = (const GroupElement & a);
	
	
	GroupElement & operator *= (const GroupElement & a);
	friend GroupElement operator * 
		(const GroupElement & a, const GroupElement & b);
		
		
	friend INT	operator == (
		const GroupElement & a, const GroupElement & b)
	{	return a.getElement()==b.getElement();	}
	friend INT	operator != (
		const GroupElement & a, const GroupElement & b)
	{	return a.getElement()!=b.getElement();	}

	
	GroupElement invert();
	void	setNeutral()
	{	element = 0;	}
	
	INT	isNeutral()
	{	return element==0;	}

	Group * getGroup() const
	{	return group;	}
	
	void	setGroup(Group *g)
	{	group = g;	}
	
	INT	getElement() const
	{	return element;	}

	void	setElement(INT i)
	{	element = i;	}


	friend ostream& operator << (ostream& s, GroupElement& g);
	friend istream& operator >> (istream& s, GroupElement& g);


private:
INT	element;
Group	*group;
};

#endif
