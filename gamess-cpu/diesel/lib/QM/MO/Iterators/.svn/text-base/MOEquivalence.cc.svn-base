//***********************************************************************
//
//	Name:			MOEquivalence.cc
//
//	Description:	set of MOs to be handled equivalently
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.04.1998
//
//
//
//
//
//***********************************************************************

#include "MOEquivalence.h"
#include "MOEquivalenceProjected.h"

#include "../../Configuration/Configuration.h"
#include "../../Configuration/ConfigurationSet.h"

#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>

using namespace std;

MOEquivalence::MOEquivalence(String equivString)
{
	n = 0;
	equiv = NULL;

	scanString(equivString);	
}


MOEquivalence::MOEquivalence(const EnergyType *e, INT maxMO)
{
INT* f = new INT[maxMO];
char	buf[100];
String	s;
	memset(f, 0, maxMO*sizeof(INT));
	
	for ( MOType i=0 ; i<maxMO ; i++ )
	{
	INT	first = 1;
		if ( f[i] )
			continue;
			
		for ( INT j=i+1 ; j<maxMO ; j++ )
		{
			if ( f[j] )
				continue;
//			cout << i+1 << " " << j+1 << " " << fabs(e[i]-e[j]) << endl;
			if ( fabs(e[i]-e[j])<degenTolerance )
			{
				f[i] = f[j] = 1;
				if ( first )
				{
					first = 0;
					sprintf(buf, "%d", i+1);
					s += String(buf);
				}
				sprintf(buf, "=%d", j+1);
				s += String(buf);
			}
		}
		if ( !first )
			s += String(" ");
		
	}
	scanString(s);
	delete[] f;
}


MOEquivalence::~MOEquivalence()
{
	if ( equiv )
	{
		for ( INT i=0 ; i<n ; i++ )
			delete equiv[i].mo;
		delete[] equiv;
	}
}


MOEquivalence::MOEquivalence(const MOEquivalence &m)
{
	n = m.n;
	equiv = new TEquiv[n];
	for ( INT i=0 ; i<n ; i++ )
	{
		equiv[i].n = m.equiv[i].n;
		equiv[i].mo = new MOType[equiv[i].n];
		memcpy(equiv[i].mo, m.equiv[i].mo, equiv[i].n*sizeof(MOType));
	}
}


void	MOEquivalence::scanString(String equivString)
{
const INT	maxn = 1000;
String	res[maxn];
	equivString.del(Regex(" *{ *"));
	equivString.del(Regex(" *} *"));
	n = split(equivString, res, maxn, Regex(" +"));

	equiv = new TEquiv[n];
	for ( INT i=0 ; i<n ; i++ )
	{
	String	res2[maxn];
		equiv[i].n = split(res[i], res2, maxn, String("="));
		equiv[i].mo = new MOType[equiv[i].n];
		for ( INT j=0 ; j<equiv[i].n ; j++ )
		{
		int h;
			sscanf(res2[j].chars(),"%d" , &h);
			equiv[i].mo[j] = h;
		}
	}
}



String	MOEquivalence::getEquivalence() const
{
String	s = " ";
char	buf[100];

	for ( INT i=0 ; i<n ; i++ )
	{
		for ( INT j=0 ; j<equiv[i].n ; j++ )
		{
			sprintf(buf, "%d=", equiv[i].mo[j]);
			s += String(buf);
		}
		s.del(((INT) s.length()-1), 1);
		s += " ";
	}
	
	return s;
}

void	MOEquivalence::symmetrize(ConfigurationSet &confSet, const MRMOs *mrmos)
{
ConfigurationSet set(confSet);

Pix	i = set.first();
	while(i)
	{
	MOEquivalenceProjected	projected(*this, set(i), mrmos);
	
		projected.generate(confSet);
	
		set.next(i);
	}

}




ostream & operator << (ostream &s, const MOEquivalence &r)
{
	s << r.getEquivalence() ;
	return s;
}

const EnergyType MOEquivalence::degenTolerance = 1e-2;
