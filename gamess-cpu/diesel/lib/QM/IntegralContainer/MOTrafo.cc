//***********************************************************************
//
//	Name:			MOTrafo.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.05.1998
//
//
//
//
//
//***********************************************************************

#include "MOTrafo.h"

#include "../../Math/MatrixVector/Matrix.h"
#include "../../Math/MatrixVector/Vector.h"
#include "../../Container/GrepAwk.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <map>

using namespace std;

MOTrafo::MOTrafo(const MOIrReps &moirreps, istream &s) :
	MOIrReps(moirreps)
{
	coefs = new MOCoefType * [IrReps];
	occNum = new MOCoefType * [IrReps];

const INT maxchars = 200;
char	buf[maxchars];


	s.getline(buf, maxchars);
	for ( INT i=0 ; i<IrReps ; i++ )
	{
		coefs[i] = new MOCoefType [inIrRep[i]*inIrRep[i]];
		occNum[i] = new MOCoefType [inIrRep[i]];
		memset(occNum[i], 0, inIrRep[i]*sizeof(MOCoefType));
		
	MOCoefType	*p = coefs[i];
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
		{
			s.getline(buf, maxchars);
			
			for ( INT k=0 ; k<inIrRep[i] ; k++ )
			{
				s.read(buf, 18);
				buf[18] = 0;
				*p++ = atof(buf);
				if ( k%4==3 || k==inIrRep[i]-1 )
					s.getline(buf, maxchars);
			}
		}
	}
}

MOTrafo::MOTrafo(istream &s) :
	MOIrReps()
{
GrepAwk	ga(s);
map<INT,INT>	M;
	while ( 1 )
	{
		if ( ga.grep("ORBITAL") )
		{
			M[atoi(ga.getWord(3).chars())] = atoi(ga.getWord(4).chars());
			++maxMO;
		}
		
		if ( ga.illegal() )
			break;
		ga++;
	}
	IrReps = M.size();
	inIrRep = new INT [IrReps];
	for ( map<INT,INT>::const_iterator i=M.begin() ; i!=M.end() ; ++i )
		inIrRep[i->first-1] = i->second;
	init();
	coefs = new MOCoefType * [IrReps];
	occNum = new MOCoefType * [IrReps];


	ga.head();
	ga++;
	
	for ( INT i=0 ; i<IrReps ; i++ )
	{
		coefs[i] = new MOCoefType [inIrRep[i]*inIrRep[i]];
		occNum[i] = new MOCoefType [inIrRep[i]];
		memset(occNum[i], 0, inIrRep[i]*sizeof(MOCoefType));
		
	MOCoefType	*p = coefs[i];
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
		{
			ga++;
			for ( INT k=0 ; k<inIrRep[i] ; k++ )
			{
				*p++ = atof(ga.getLine().at(k%4*18, 18).chars());
				if ( k%4==3 || k==inIrRep[i]-1 )
					ga++;
			}
		}
	}
	ga++;
	for ( INT i=0 ; i<IrReps ; i++ )
	{
	MOCoefType	*p = occNum[i];
		for ( INT k=0 ; k<inIrRep[i] ; k++ )
		{
			*p++ = atof(ga.getLine().at(k%4*18, 18).chars());
			if ( k%4==3 || k==inIrRep[i]-1 )
				ga++;
		}
	}
}



MOTrafo::~MOTrafo()
{
	if ( coefs )
	{
		for ( INT i=0 ; i<IrReps ; i++ )
			if ( coefs[i] )
				delete coefs[i];
		delete coefs;
	}

	if ( occNum )
	{
		for ( INT i=0 ; i<IrReps ; i++ )
			if ( occNum[i] )
				delete occNum[i];
		delete occNum;
	}
}


MOTrafo::operator Matrix<MOTrafo::MOCoefType> () const
{
MOCoefType* pp = new MOCoefType[maxMO*maxMO];

	unpack(pp);
	
Matrix<MOCoefType>	prod(maxMO, maxMO, pp);
	return prod;

delete[] pp;
}

void	MOTrafo::transform(const Matrix<MOCoefType> &mat)
{
MOCoefType* pp = new MOCoefType[maxMO*maxMO];

	unpack(pp);
	
Matrix<MOCoefType>	prod(maxMO, maxMO, pp);
	prod = mat*prod;		

	pack(prod.getP());	

delete[] pp;
}

void	MOTrafo::setOccupationNumbers(const Vector<MOCoefType> &v)
{
INT	k = 0;
	for ( IrRep i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
			occNum[i][j] = v[k++];
}

Vector<MOTrafo::MOCoefType>	MOTrafo::getOccupationNumbers() const
{
Vector<MOCoefType>	v(getMaxMO());
INT	k = 0;
	for ( IrRep i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
			v[k++] = occNum[i][j];
	return v;
}



void	MOTrafo::unpack(MOCoefType *pp) const
{
	memset(pp, 0, maxMO*maxMO*sizeof(MOCoefType));
	for ( INT i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
			for ( INT k=0 ; k<inIrRep[i] ; k++ )
				pp[
					(getStartMO(i)+j-1)*maxMO +
					(getStartMO(i)+k-1)
					] = coefs[i][j*inIrRep[i] + k];
}


void	MOTrafo::pack(const MOCoefType *pp)
{
	for ( INT i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
			for ( INT k=0 ; k<inIrRep[i] ; k++ )
				coefs[i][j*inIrRep[i] + k] = pp[
					(getStartMO(i)+j-1)*maxMO +
					(getStartMO(i)+k-1)
					];
}

ostream & operator << (ostream &os, const MOTrafo &t)
{
MOTrafo::MOCoefType* pp = new MOTrafo::MOCoefType[t.maxMO*t.maxMO];

	t.unpack(pp);
	
	for ( INT i=0 ; i<t.maxMO ; i++ )
	{
		for ( INT j=0 ; j<t.maxMO ; j++ )
			os << setw(12) << pp[j*t.maxMO + i] << " ";
		os << endl;
	}

	delete[] pp;

	return os;
}


ostream & MOTrafo::writeToStream(ostream &s) const
{
	s << "* Natural Orbitals" << endl;
	s.precision(11);
	s.setf(ios::scientific);
	for ( IrRep i=0 ; i<IrReps ; i++ )
	{
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
		{
			s << "* ORBITAL" << setw(5) << i+1 << setw(5) << j+1 << endl;
			for ( INT k=0 ; k<inIrRep[i] ; k++ )
			{
				s << setw(18) << (*this)(i, j, k);
				if ( k%4==3 )
					s << endl;
			}
			if ( inIrRep[i]%4 )
				s << endl;
		}	
	}
	s << "* OCCUPATION NUMBERS" << endl;
	for ( IrRep i=0 ; i<IrReps ; i++ )
	{
		for ( INT j=0 ; j<inIrRep[i] ; j++ )
		{
			s << setw(18) << occNum[i][j];
			if ( j%4==3 )
				s << endl;
		}
		if ( inIrRep[i]%4 )
			s << endl;
	}
	return s;
}

