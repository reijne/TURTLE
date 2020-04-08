//***********************************************************************
//
//	Name:			Group.cc
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
#include "Group.h"

#define EPS 1E-6

using namespace std;

//*********************************************************************
//
// Implementation Group
//


Group::Group()
{
	elements = 0;
}

Group::Group(INT *prodTab, INT _elements)
{
	elements = _elements;
	if ( elements )
	{	ProdTab = new INT[elements*elements];
		memcpy(ProdTab, prodTab, elements*elements*sizeof(INT));
		invTab = new INT[elements];
		initInvTab();
		initClasses();
		initIrredMult();
		initCharTab();
	}
}


Group::~Group()
{
	if ( elements )
	{	elements = 0;
		delete ProdTab;
		delete invTab;
		delete nPerClass;
		delete[] CharTab;
		for ( INT i=0 ; i<nClasses ; i++ )
			delete Classes[i];
		delete Classes;
		nClasses = 0;
	}
}

void Group::setGroup(INT *prodTab, INT _elements)
{
	if ( elements )
	{	elements = 0;
		delete ProdTab;
		delete invTab;
		delete nPerClass;
		delete[] CharTab;
//		for ( INT i=0 ; i<nClasses ; i++ )
//			delete Classes[i];
		delete Classes;
		nClasses = 0;
	}
	elements = _elements;
	if ( elements )
	{	ProdTab = new INT[elements*elements];
		memcpy(ProdTab, prodTab, elements*elements*sizeof(INT));
		invTab = new INT[elements];
		initInvTab();
		initClasses();
		initIrredMult();
		initCharTab();
	}
}


void	Group::setNumberOfElements(INT n)
{
	if ( n==elements )
		return;
	if ( elements )
	{	elements = 0;
		delete ProdTab;
		delete invTab;
		delete nPerClass;
		delete[] CharTab;
//		for ( INT i=0 ; i<nClasses ; i++ )
//			delete Classes[i];
		delete Classes;
		nClasses = 0;
	}
	elements = n;
	if ( elements )
		ProdTab = new INT[elements*elements];
}



void	Group::initInvTab()
{
	for ( INT i=0 ; i<elements ; i++ )
	{	INT n = i*elements;
		for ( INT j=0 ; j<=elements ; j++ )
			if ( !ProdTab[n+j] )
			{	invTab[i] = j;
				break;
			}
	}
}

void	Group::initClasses()
{
        INT* h = new INT[elements];
	for ( INT i=0 ; i<elements ; i++ )
		h[i] = -1;

	nClasses = 0;
INT	k;
	for ( INT i=0 ; i<elements ; i++ )
	{	if ( h[i]!= -1 )
			continue;
		for ( INT j=0 ; j<elements ; j++ )
		{	// C = A^-1 * B * A
			k = ProdTab[invTab[j]*elements + ProdTab[i*elements + j]];
			h[k] = nClasses;
		}
		nClasses++;
	}
//	for ( INT i=0 ; i<elements ; i++ )
//		cout << i << " " << h[i] << endl;
		
		
	nPerClass = new INT[nClasses];
	Classes = new INT * [nClasses];
	for ( INT i=0 ; i<nClasses ; i++ )
	{	nPerClass[i] = 0;
		for ( INT j=0 ; j<elements ; j++ )
			if ( h[j] == i )
				nPerClass[i]++;
//		cout << "nPerClass: " << nPerClass[i] << endl;
		Classes[i] = new INT[nPerClass[i]];
		k = 0;
		for ( INT j=0 ; j<elements ; j++ )
			if ( h[j] == i )
				Classes[i][k++] = j;	
	}
	delete h;
}


void	Group::initIrredMult()
{
INT	*IrredMult = new INT[nClasses];

	CharTab = new Complex[nClasses*nClasses];
	
	memset(CharTab, 0, nClasses*nClasses*sizeof(Complex));
	
	calcIrredMult(nClasses-1, IrredMult+nClasses-1, elements - nClasses, 3);

	for ( INT i=0 ; i<nClasses ; i++ )
	{	cout << "Irred: " << i << ", Mult " << IrredMult[i] << endl;
		CharTab[i*nClasses] = IrredMult[i];
		CharTab[i] = 1;
	}

	delete IrredMult;
		
}

void	Group::initCharTab()
{
INT	n = 0;

Complex* a = new Complex[nClasses];
INT	i;

	i = 0;
	for ( INT j=0 ; j<nClasses ; j++ )
		a[j] = 1;
		
	for ( ; i>=0 ; )
	{
/*		for ( INT j=0 ; j<nClasses ; j++ )
			cout << a[j];
		cout << endl;
*/		
		if ( checkOrtho(a, n) )
		{	memcpy(&CharTab[n*nClasses], a, nClasses*sizeof(Complex));
			n++;
		}
		
		
		i = nClasses-1;
		while ( i+1 )
		{	if ( a[i] == -1.0 )
			{	a[i] = 1;
				i--;
				continue;
			}
			a[i] -= 2;
			break;
		} 
	}



/*
	for ( INT i=1 ; i<nClasses ; i++ )
	{
	
	
	
		INT Mult = (INT) real(CharTab[i*nClasses]);
		if ( oldMult!=Mult )
		{	oldMult = Mult;
			n = calcBetraege(Mult, Betraege);
		}
		cout << endl;
	}
*/


	for ( i=0 ; i<nClasses ; i++ )
	{	for ( INT j=0 ; j<nClasses ; j++ )
			cout << CharTab[i*nClasses+j] << "\t";
		cout << endl;
	}
delete a;
}


INT	Group::checkOrtho(Complex *a, INT n)
{
INT	ortho = 1;
	for ( INT j=0 ; j<n ; j++ )
		if ( norm(calcProd(&CharTab[j*nClasses], a))>EPS )
		{	ortho = 0;
			break;
		}
	return ortho;
}


void	Group::calcBetraege(INT Mult, INT *Betraege)
{
INT	n = 0;
const INT maxN = 3;
INT* ii = new INT[nClasses];
INT	i = 1;
INT	rest = elements - Mult*Mult;
	ii[0] = Mult;
	ii[i]=0;
	for ( ; ; )
	{	if ( ii[i]>maxN )
		{	i--;
			rest += ii[i]*ii[i]*nPerClass[i];
			ii[i]++;
			continue;
		}
		rest -= ii[i]*ii[i]*nPerClass[i];
//		cout << "i= " << i << ", ii[i]= " << ii[i] << ", rest= " << rest << endl;
		if ( rest>=0 )
		{	if ( i<nClasses-1 )
			{	i++;
				ii[i]=0;
				continue;
			}
			else
			if ( rest==0 )
			{	memcpy(&Betraege[n], ii, nClasses*sizeof(INT));
				for ( i=0 ; i<nClasses ; i++ )
					cout << ii[i] << "\t";
				cout << endl;
				n++;
				if ( n==nClasses )
					break;
			}
		}
		rest += ii[i]*ii[i]*nPerClass[i];
		if ( ii[i]>=maxN )
		{	i--;
			if ( i<=0 )
				break;
			rest += ii[i]*ii[i]*nPerClass[i];
		}
		
		ii[i]++;
	}
delete ii;
}

INT	Group::calcBetrag(INT *Betraege)
{
INT	Betrag = 0;

	for ( INT i=0 ; i<nClasses ; i++ )
		Betrag += Betraege[i]*Betraege[i];
		
	return Betrag;
}

Complex	Group::calcProd(Complex *p1, Complex *p2)
{
Complex	result = 0;

	for ( INT i=0 ; i<nClasses ; i++ )
		result += conj(p1[i])*p2[i]*(double)nPerClass[i];
		
	return result;
}


INT	Group::calcIrredMult(INT n, INT *p, INT rest, INT max)
{
INT j;
	if ( rest<0 )
		return 0;

	for ( j=2 ; j<=max ; j++ )
		if ( rest-j*j+1<0 )
			break;
	do
	{	j--;
		*p = j;
		if ( n )
		{	if ( calcIrredMult(n-1, p-1, rest-j*j+1, j) )
				return 1;
		}
		else
			return (rest-j*j+1==0);
		
	} while ( j>1 );
	return 0;
}

ostream& operator<<(ostream& s, Group & a)
{
INT	elements = a.getNumberOfElements();
INT	i, j;

	s << "    ";
	for ( i=0 ; i<elements ; i++ )
		s << i << " ";
	s << "\n--+";
	for ( i=0 ; i<elements ; i++ )
		s << "--";
	s << "\n";
	for ( i=0 ; i<elements ; i++ )
	{	s << i << " | ";
		for ( j=0 ; j<elements ; j++ )
			s << a[i*elements + j] << " ";
		s << "\n";
	}
	return s;	
}

/*
istream& operator>>(istream& s, Group & a)
{
}
*/
