//***********************************************************************
//
//	Name:			SelIterator.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			28.10.1998
//
//
//
//
//
//***********************************************************************




#include "SelIterator.h"


#include "../MRTree/Diag/NExternalsDiag.h"
#include "../MRTree/Diag/InternalConfsDiag.h"
#include "../MRTree/Diag/TupelStructureDiag.h"
#include "../MRTree/Diag/extMOsDiag.h"
#include "../MRTree/Container/ContainerIterator.h"

#include "../MO/Iterators/MOAccess.h"

#include "../../Container/AVLSet.h"
#include "../Configuration/ConfigurationStart.h"


SelIterator::SelIterator(
	const NExternalsDiag & Psi,
	const AVLSet<ConfigurationStart<MOType> > *internal,
	MOAccess	* const *moAccess,
	INT maxExc) :
	Psi(Psi),
	internal(internal),
	moAccess(moAccess),
	maxExc(maxExc)
{
	end = 0;
	confStart = 0;
	SAFStart = 0;
	SAFInc = 0;
	
	nExt = 0;
	p0 = NULL;
	p1 = NULL;
	p2 = NULL;
	next(1);
}

SelIterator::~SelIterator()
{
}



void	SelIterator::next(INT first)
{
start:
	switch ( nExt ) {
	case 0:
		break;

	case 1:
		if ( iter3<p3->getNumberOfElements() )
		{
			SAFStart = lSAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0]);
			confStart = lconfStart + moAccess[nExt]->getConfOffset(mo[0]);
			iter3++;
			mo++;
			return;
		}
		break;

	case 2:
		if ( iter3<p3->getNumberOfElements() )
		{
			if ( p3->getNumberOfOpenMOs()==2 )
			{
				SAFStart = lSAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[1]);
				confStart = lconfStart + moAccess[nExt]->getConfOffset(mo[0], mo[1]);
				iter3++;
				mo += 2;
				return;
			}
			else
			{
				SAFStart = lSAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[0]);
				confStart = lconfStart + moAccess[nExt]->getConfOffset(mo[0], mo[0]);
				iter3++;
				mo++;
				return;
			}
		}
		break;
	}


INT	prepairN = first;

	if ( !p0 )
	{
		p0 = &Psi;
		iter0 = p0->first();
		nExt = 0;
	}
	while ( !p0->isLast(iter0) )
	{
		if ( !p1 )
		{
			p1 = (*p0)[iter0];
			iter1 = p1->first();
		}
		while ( !p1->isLast(iter1) )
		{
			if ( !p2 )
			{
				p2 = (*p1)[iter1];
				iter2 = p2->first();
			}
			while ( !p2->isLast(iter2) )
			{
				p3 = (*p2)[iter2];
				if ( !p3 || !p3->getNumberOfElements() )
				{
					p2->next(iter2);
					continue;
				}
				iter3 = 0;
				prepair(prepairN);
				p2->next(iter2);
				if ( nExt )
					goto start;
				else
					return;
			}
			p1->next(iter1);
			p2 = NULL;
			prepairN = 1;
		}
		p0->next(iter0);
		p1 = NULL;
		nExt++;
	}
	end = 1;
}



void	SelIterator::prepair(INT prepairN)
{
	if ( prepairN )
	{
	ConfigurationStart<MOType>	conf1(*p2);
	Pix	pix = internal[nExt].seek(conf1);
		if ( !pix )
		{
			cout << "wrong Psi wave function" << endl;
			cout << ((Configuration<MOType>) conf1) << endl;
			exit(1);
		}
		lSAFStart = internal[nExt](pix).getSAFStart();
		lconfStart = internal[nExt](pix).getConfStart();
		internalOpen = internal[nExt](pix).getNumberOfOpenShells();
	}

	switch ( nExt ) {
	case 0:
		SAFStart = lSAFStart;
		confStart = lconfStart;
		SAFInc = Psi.getNumberOfSpinAdaptedFunctions(internalOpen);
		mo = NULL;
		break;

	case 1:
		SAFInc = Psi.getNumberOfSpinAdaptedFunctions(internalOpen+1);
		mo = (*((extMOsDiag *) p3))[0];
		break;

	case 2:
		if ( p3->getNumberOfOpenMOs()==2 )
			SAFInc = Psi.getNumberOfSpinAdaptedFunctions(internalOpen+2);
		else
			SAFInc = Psi.getNumberOfSpinAdaptedFunctions(internalOpen);
		mo = (*((extMOsDiag *) p3))[0];
		break;
	}
}


/*
void	SelIterator::next()
{

	for ( iter0 = Psi.first() ; !Psi.isLast(iter0) ; Psi.next(iter0) , nExt++ )
	{
		p1 = Psi[iter0];

		for ( iter1 = p1->first() ; !p1->isLast(iter1) ; p1->next(iter1) )
		{
			p2 = (*p1)[iter1];

			if ( !p2->getNumberOfElements() )
				continue;

		ConfigurationStart<MOType>	conf1(*p2);
		Pix	pix = internal[nExt].seek(conf1);
			if ( !pix )
			{
				cout << "wrong Psi wave function" << endl;
				cout << ((Configuration<MOType>) conf1) << endl;
				exit(1);
			}

		INT	SAFStart = internal[nExt](pix).getSAFStart();
		INT	internalOpen = internal[nExt](pix).getNumberOfOpenShells();
//				cout << "::::::::::::::::  SAFStart: " << SAFStart << endl;

			for ( iter2 = p2->first() ; !p2->isLast(iter2) ; p2->next(iter2) )
			{
				p3 = (*p2)[iter2];
				if ( !p3->getNumberOfElements() )
					continue;


				switch ( nExt ) {
				case 0:
					{
//							cout << SAFStart << "    0" << endl;

					INT	safs = Psi.getNumberOfSpinAdaptedFunctions(internalOpen);
//						memset(&v[SAFStart], 0, safs*sizeof(VectorType));
					}

					break;

				case 1:
					{
						mo = (*((extMOsDiag *) p3))[0];
					INT	safs = Psi.getNumberOfSpinAdaptedFunctions(internalOpen+1);
						for ( iter3 = 0 ; iter3<p3->getNumberOfElements() ; iter3++ , mo++ )
						{	
//								cout << SAFStart  + moAccess[nExt]->getSAFOffset(internalOpen, mo[0])<< "    1" << endl;
//							memset(&v[SAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0])], 0, safs*sizeof(VectorType));
						}
					}
					break;

				case 2:
					{
						mo = (*((extMOsDiag *) p3))[0];
						if ( p3->getNumberOfOpenMOs()==2 )
						{
						INT	safs = Psi.getNumberOfSpinAdaptedFunctions(internalOpen+2);
							for ( iter3 = 0 ; iter3<p3->getNumberOfElements() ; iter3++ , mo+=2 )
							{	
//									cout << SAFStart  + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[1])<< "    2" << endl;
//								memset(&v[SAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[1])], 0, safs*sizeof(VectorType));
							}
						}
						else
						{
						INT	safs = Psi.getNumberOfSpinAdaptedFunctions(internalOpen);
							for ( iter3 = 0 ; iter3<p3->getNumberOfElements() ; iter3++ , mo++ )
							{	
//									cout << SAFStart  + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[0])<< "    2" << endl;
//								memset(&v[SAFStart + moAccess[nExt]->getSAFOffset(internalOpen, mo[0], mo[0])], 0, safs*sizeof(VectorType));
							}
						}
					}

					break;

				}
			}

		}
	}
}
*/
