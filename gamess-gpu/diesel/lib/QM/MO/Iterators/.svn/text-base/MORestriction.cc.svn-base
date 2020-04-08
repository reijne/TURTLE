//***********************************************************************
//
//	Name:			MORestriction.cc
//
//	Description:	restrictions for MO occupation patterns
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.03.1998
//
//
//
//
//
//***********************************************************************

#include "MORestriction.h"
#include "../MOIrReps.h"

#include "../../Configuration/Configuration.h"

#include <ctype.h>
#include <sstream>
#include <stdio.h>

//#include <mask>

using std::istringstream;

char *MORestriction::Operators[] = {
		"<", "<=", "=", ">", ">=", "!=" 
	};


MORestriction::MORestriction(INT _maxMO, INT _nRestrictions)
{
	maxMO = _maxMO;
	nRestrictions = _nRestrictions;
	restrictions = new TRestriction[nRestrictions];
}

MORestriction::MORestriction(
	const Configuration<MOType> &conf, 
	const MOIrReps &moirreps,
	INT activeVirtual)
{
INT	nIrreps = moirreps.getNumberOfIrReps();
INT* highest = new INT[nIrreps];
	for ( IrRep i=0 ; i<nIrreps ; i++ )
		highest[i] = moirreps.getStartMO(i) - 1;
	
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		highest[moirreps.getIrRep(conf.getOpenShell(i))] =
			conf.getOpenShell(i);
	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
	{
		if ( highest[moirreps.getIrRep(conf.getClosedShell(i))] < 
					conf.getClosedShell(i) )
			highest[moirreps.getIrRep(conf.getClosedShell(i))] =
				conf.getClosedShell(i);
	}
	
	nRestrictions = nIrreps;			
	restrictions = new TRestriction[nIrreps];
	maxMO = moirreps.getMaxMO();

	for ( IrRep i=0 ; i<nIrreps ; i++ )
	{
		restrictions[i].occupation = 0;
		restrictions[i].op = 2;	// "="
		
		for ( MOType j = highest[i]+activeVirtual ; j<=moirreps.getEndMO(i) ; j++ )
			restrictions[i].mask.set(j-1);
	}
	delete[] highest;
		

}

void	MORestriction::setMaxMO(INT _maxMO)
{
	maxMO = _maxMO;
}


MORestriction::MORestriction(
	String restString)
{
	maxMO = 0;
	
String	s;

	nRestrictions = 0;
	{
	INT j=0;
		while ( (j=restString.index(Regex("[<>!=]=?"), j+1)+1)>0 )
			nRestrictions++;
	}
	restrictions = new TRestriction[nRestrictions];


	restString += " ";
INT	i = 0;
INT	n = 0;
	while ( i<(INT) restString.length() )
	{
	String	s(String(restString.from(i)).before(' '));
		if ( s.length()>=3 )
		{
		INT j = s.length()-1;
			while ( j>=0 && isdigit(s[j]) )
				j--;
			while ( j>=0 && !isdigit(s[j]) )
				j--;
			j++;
		
			getOperatorOccupation(s.from(j).chars(), 
				restrictions[n].op, restrictions[n].occupation);
		String	range(substRange(s.before(j)));
		istringstream	str(range.chars(), istringstream::in | istringstream::out);
		MOType	mo;
			while ( str>>mo )
			{
				restrictions[n].mask.set(mo-1);
				if ( mo>maxMO )
					maxMO = mo;
			}
			n++;
		}
		i += s.length()+1;
	}
}


MORestriction::~MORestriction()
{
	delete[] restrictions;
}

MORestriction::MORestriction(const MORestriction &m)
{
	maxMO = m.maxMO;
	nRestrictions = m.nRestrictions;
	restrictions = new TRestriction[nRestrictions];
	for ( INT i=0 ; i<nRestrictions ; i++ )
	{
		restrictions[i].mask = m.restrictions[i].mask;
		restrictions[i].occupation = m.restrictions[i].occupation;
		restrictions[i].op = m.restrictions[i].op;
	}
}


INT	MORestriction::check(
	const Configuration<MOType> &conf) const
{
BitSet	doubleMOs, singleMOs;

	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		singleMOs.set(conf.getOpenShell(i)-1);

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		doubleMOs.set(conf.getClosedShell(i)-1);


INT	occupation;

	for ( INT i=0 ; i<nRestrictions ; i++ )
	{
		occupation = 
			2*(doubleMOs & restrictions[i].mask).count() +
			(singleMOs & restrictions[i].mask).count();

		switch ( restrictions[i].op ) {
		case 0:
			if ( !(occupation<restrictions[i].occupation) )
				return 0;
			break;
		
		case 1:
			if ( !(occupation<=restrictions[i].occupation) )
				return 0;
			break;
		
		case 2:
			if ( !(occupation==restrictions[i].occupation) )
				return 0;
			break;
		
		case 3:
			if ( !(occupation>restrictions[i].occupation) )
				return 0;
			break;
		
		case 4:
			if ( !(occupation>=restrictions[i].occupation) )
				return 0;
			break;
		
		case 5:
			if ( !(occupation!=restrictions[i].occupation) )
				return 0;
			break;
		}
	}
	return 1;

}



void	MORestriction::getOperatorOccupation(
	const char *p,
	INT	&op, INT &occupation) const
{
INT	l = 1;
	if ( p[1]=='=' )
		l++;
	
	op = occupation = -1;
	for ( INT i=0 ; i<nOperators ; i++ )
		if ( !strncmp(Operators[i], p, l) )
		{
			op = i;
			occupation = atoi(p+l);
			return;
		}
	return;
}

String	MORestriction::substRange(String s, char sep) const
{
INT	i = 0;
INT	j;
	while ( (j=s.index('-', i))>=0 )
	{
String	result;
	INT	from = j-1;
		while ( from>=0 && isdigit(s[from]) )
			from--;
		from++;
	INT	to = j+1;
		while ( to<(INT) s.length() && isdigit(s[to]) )
			to++;
		to--;
		
		for ( INT k=atoi(s.at(from, j-from).chars()) ; 
			k<=atoi(s.at(j+1, to-j).chars()) ; k++ )
		{
		char	buf[10];
			sprintf(buf, "%d", k);
			result +=  String(buf) + sep;
		}
		s = s.before(from) + result + s.after(to);
		i = j+1;
	}
	s.gsub(",", " ");
	return s;
}




String	MORestriction::getRestriction() const
{
String	s = " ";
char	buf[100];

	for ( INT i=0 ; i<nRestrictions ; i++ )
	{
	INT	hit = 0;
	INT start = -1;
		for ( MOType mo=0 ; mo<maxMO ; mo++ )
		{
		INT	miss = restrictions[i].mask.test(mo);
				
			if ( miss )
			{
				if ( start<0 )
					start = mo;
				hit = 1;
			}
			if ( !miss || mo==maxMO-1 )
			{
				if ( start>=0 )
				{
					if ( mo-2>start )
					{
						sprintf(buf, "%d-%d,", start+1, mo-1+1+(mo==maxMO-1));
						s += String(buf);
					}
					else
					if ( mo-1>start )
					{
						sprintf(buf, "%d,%d,", start+1, mo-1+1+(mo==maxMO-1));
						s += String(buf);
					}
					else
					{
						sprintf(buf, "%d,", start+1);
						s += String(buf);
					}
				}
				start = -1;
			}
		}
		if ( hit )
		{
			s.del(((INT) s.length()-1), 1);
			sprintf(buf, "%s%d ", Operators[restrictions[i].op],
				restrictions[i].occupation);
			s += String(buf);
		}
		
	}
	return s;
}



INT	MORestriction::isIllegal() const
{
//	for ( MOType mo=0 ; mo<maxMO ; mo++ )
//		if ( mask[0].test(mo) && mask[1].test(mo) && mask[2].test(mo) )
//			return 1;
	return 0;
}

/*MORestriction &MORestriction::operator |= (const MORestriction &r)
{
	if ( r.maxMO>maxMO )
		maxMO = r.maxMO;
		
	mask[0] |= r.mask[0];
	mask[1] |= r.mask[1];
	mask[2] |= r.mask[2];
	return *this;
}
*/


ostream & operator << (ostream &s, const MORestriction &r)
{
	s << r.getRestriction() ;
	return s;
}


