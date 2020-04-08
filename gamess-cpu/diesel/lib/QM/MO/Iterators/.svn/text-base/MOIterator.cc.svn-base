//***********************************************************************
//
//	Name:			MOIterator.cc
//
//	Description:	MOs classified by internal/external status and
//					irrep to be used as creators or annihilators on
//					internal configuration rests
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			20.09.1998
//
//
//
//
//
//***********************************************************************

#include "MOIterator.h"

#include <iostream>

using namespace std;

MOIterator::MOIterator(
	INT nTupel, const MRMOs *mrmos, INT isExternal) :
	nTupel(nTupel), mrmos(mrmos), isExternal(isExternal), 
	specificIrrep(0)
{
	init();
}



MOIterator::MOIterator(
	INT nTupel, const MRMOs *mrmos, INT isExternal, IrRep productIrRep) :
	nTupel(nTupel), mrmos(mrmos), isExternal(isExternal), 
	productIrRep(productIrRep),
	specificIrrep(1)
{
	init();
}





MOIterator::~MOIterator()
{
	if ( irrep )
		delete[] irrep;
	if ( mos )
		delete[] mos;
}


MOIterator::operator Configuration<MOType>() const
{
Configuration<MOType>	conf;
	for ( INT i=0 ; i<nTupel ; i++ )
		conf.create(mos[i]);

	return conf;
}


void	MOIterator::init()
{
	irrep = NULL;
	mos = NULL;
	if ( nTupel )
	{
		changed = 1;
		irrep = new IrRep[nTupel];
		irrep[0] = -1;
		lIrRep = 0;
		mos = new MOType[nTupel];
		do {
			if ( (end = !irrepIter()) )
				return;

			mos[0] = mrmos->getIrRepIntExtStart(irrep[0], isExternal) - 1;
			iMOs = 0;
			if ( moIter() )
				return;
		} while ( 1 );
	}
	else
	{
		end = 0;
	}
}


INT	MOIterator::irrepIter()
{
	irrep[lIrRep]++;
	//	loops over symmetries
	for ( ; ; )
	{
		if ( irrep[lIrRep] >= mrmos->getNumberOfIrReps() )
		{
			lIrRep--;
			if ( lIrRep<0 )
				return 0;
			else
			{
				irrep[lIrRep]++;
				continue;
			}
		}
		if ( lIrRep<nTupel-1 )
		{
			lIrRep++;
			irrep[lIrRep] = irrep[lIrRep-1];
			continue;
		}
		if ( specificIrrep )
		{
			//	check matching symmetry
		IrRep prod = 0;
			for ( INT i=0 ; i<nTupel ; i++ )
			{
//				cout << irrep[i] << " ";
				prod = mrmos->getProd(prod, irrep[i]);
			}

//			cout << "prod=" << prod << endl;
			if ( prod != productIrRep )
			{
				irrep[lIrRep]++;
				continue;
			}
		}
		return 1;
	}	//	loops over symmetries
}

void	MOIterator::skip(INT n)
{
	for ( INT i=0 ; i<n ; i++ )
		if ( !moIter() )
			break;
}

INT	MOIterator::moIter()
{

	mos[iMOs]++;

	// loop over MOs
	for ( ; ; )
	{
		if ( mos[iMOs] > getMOEnd(iMOs) )
		{
			iMOs--;
			if ( iMOs<0 )
				return 0;
			else
			{
				mos[iMOs]++;
				continue;
			}
		}
		
		if ( iMOs<nTupel-1 )
		{
			iMOs++;
			mos[iMOs] = getMO(mos[iMOs-1], iMOs);
			continue;
		}
		return 1;
	}

}


void	MOIterator::next()
{
	end |= !nTupel;
	if ( end )
		return;
	changed = 0;
		

	if ( moIter() )
		return;


	// symmetry changed ==> interaction case may have changed
	changed = 1;
	do
	{
		if ( (end = !irrepIter()) )
			return;

		mos[0] = mrmos->getIrRepIntExtStart(irrep[0], isExternal) - 1;
		iMOs = 0;
	} while ( (end = !moIter()) );


}



ostream & operator << (ostream & s, const MOIterator & m)
{
	for ( INT i=0 ; i<m.nTupel ; i++ )
		s << m.irrep[i] << " ";
	s << endl;
	

	s << m.operator Configuration<MOType>() << endl;
	return s;
}
