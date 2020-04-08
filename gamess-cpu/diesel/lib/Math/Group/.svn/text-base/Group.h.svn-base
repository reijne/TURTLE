//***********************************************************************
//
//	Name:			Group.h
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

#ifndef _GROUP_H
#define _GROUP_H

#include "../../../config.h"

#include <iostream>
using std::ostream;
using std::istream;

#include "Complex.h"

class Group {
public:
	Group();
	Group(INT *prodTab, INT Elements);
	~Group();
	
	
	INT	getProdukt(const INT & a, const INT & b)
	{	return ProdTab[a*elements + b];	}
	
	INT	getInv(INT a)
	{	return	invTab[a];	}
	
	
	INT	operator[](INT n)
	{	return ProdTab[n];	}
	
	friend ostream& operator << (ostream& s, Group& g);
	friend istream& operator >> (istream& s, Group& g);

	INT	getNumberOfElements()
	{	return elements;	}
	
	INT	*getProdTab()
	{	return ProdTab;	}

protected:
	void	setGroup(INT *prodTab, INT Elements);
	void	setNumberOfElements(INT n);

	void	initInvTab();

	void	initClasses();

	void	initIrredMult();
	INT	calcIrredMult(INT n, INT *p, INT rest, INT max);

	void	initCharTab();
	void	calcBetraege(INT Mult, INT *Betraege);
	INT	calcBetrag(INT *Betraege);
	Complex	calcProd(Complex *p1, Complex *p2);
	INT	checkOrtho(Complex *a, INT n);
	

private:
INT	elements;
INT	*ProdTab;
INT	*invTab;
INT	**Classes;
INT	*nPerClass;
INT nClasses;
Complex	*CharTab;
};

#endif
