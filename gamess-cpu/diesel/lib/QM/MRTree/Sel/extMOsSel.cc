//***********************************************************************
//
//	Name:			extMOsSel.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************

#include "NExternalsSel.h"
#include "InternalConfsSel.h"
#include "TupelStructureSel.h"
#include "extMOsSel.h"

#include <stdlib.h>
#include <string>

#include <math.h>

#include "../../Configuration/DiffConf.h"


/*
extMOsSel::extMOsSel(istream &s) :
	extMOsBase(s),
	ListContainer<extEntry>(s),
	Tree<TupelStructureBase>(NULL)
{
//	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
//		operator [] (iter)->setParent(this);
}
*/


extMOsSel::extMOsSel(INT _Syms, INT _open, INT _closed, 
	TupelStructureBase *parent) :
        Tree<TupelStructureBase>(parent),
	extMOsBase(_Syms, _open, _closed)
{
	SAFInc = getParent()->getParent()->getParent()->
			getNumberOfSpinAdaptedFunctions(
			_open + getParent()->getNumberOfOpenShells());
}


extMOsSel::~extMOsSel()
{
}

extMOsSel::extMOsSel(const extMOsSel & ext) :
	extMOsBase(ext)
{
	cout << "copy constructor extMOsSel" << endl;
}

INT	extMOsSel::init(INT n, const IrRep *irrep)
{
	// initialize SAFStart
	SAFStart = n;
	SAFInc = getParent()->getParent()->getParent()->
			getNumberOfSpinAdaptedFunctions(
			open + getParent()->getNumberOfOpenShells());
	return n + length()*SAFInc;
}


Pix	extMOsSel::seek(Configuration<MOType> conf) const
{
DiffConf<MOType>	diff;

	diff.calcDiffConf(conf, *getParent());
//	cout << diff << endl;
	if ( diff.getTo().getNumberOfElectrons() )
		return NULL;
const Configuration<MOType> *rest = &diff.getFrom();

	if ( rest->getNumberOfOpenShells()!=getNumberOfOpenMOs() ||
		 rest->getNumberOfClosedShells()!=getNumberOfClosedMOs() )
		 return NULL;


ContainerIterator	iter = first();
Configuration<MOType> probe;

	while ( iter.pix )
	{
		probe.clear();
		for ( INT i=0 ; i<getNumberOfOpenMOs() ; i++ )
			probe.appendOpenShell(getOpenMO(iter.pix, i));

		for ( INT i=0 ; i<getNumberOfClosedMOs() ; i++ )
			probe.appendClosedShell(getClosedMO(iter.pix, i));

//		cout << "probe=" << probe << endl;
		if ( probe == *rest )
			return iter.pix;
		
		next(iter);
	}
	return NULL;
}


ConfigurationSAFNr<MOType>	extMOsSel::
	getConfigurationSAFNr(INT ind) const
{
ContainerIterator	iter = first();

	for ( INT i=1 ; i<ind && !isLast(iter) ; i++ )
		next(iter);

	return getConfigurationSAFNr(iter.pix);
}


ConfigurationSAFNr<MOType>	extMOsSel::getConfigurationSAFNr(Pix pix) const
{
Configuration<MOType>	conf = *getParent();


//	cout << "n= " << n << ", this = " << this << endl;
	for ( INT i=0 ; i<getNumberOfOpenMOs() ; i++ )
		conf.create(getOpenMO(pix, i));

	for ( INT i=0 ; i<getNumberOfClosedMOs() ; i++ )
	{	conf.create(getClosedMO(pix, i));
		conf.create(getClosedMO(pix, i));
	}


INT	j = 0;

Pix	ii = SLList<extEntry *>::first();
	while ( ii )
	{
		if ( ii==pix )
			break;
		SLList<extEntry *>::next(ii);
		j++;
	}


	return ConfigurationSAFNr<MOType>(
		conf, *getParent(),
		getNumberOfOpenMOs(), getNumberOfClosedMOs(),
		SAFStart + j*SAFInc, SAFInc,
		getEnergy(pix),
		getParent()->isReference(),
		getSelectionFlags(pix)
		);

}

void	extMOsSel::writeToStream(ostream & s) const
{
//	s << "Container<extEntry>::writeToStream(s)" << endl;
	Container<extEntry>::writeToStream(s);
//	s << "extMOsBase::writeToStream(s)" << endl;
	extMOsBase::writeToStream(s);
}

ostream& operator<<(ostream & s, extMOsSel & base)
{
	s << ((extMOsBase &) base);
	return s;	


/*
	for ( INT i=0 ; i<mo.getNumberOfConfs() ; i++ )
	{	for ( INT j=0 ; j<mo.getNumberOfOpenMOs() ; j++ )
			s << mo.getOpenMO(i, j) << " ";
		if ( mo.getNumberOfOpenMOs() && mo.getNumberOfClosedMOs() )
			s << "| ";
		for ( INT j=0 ; j<mo.getNumberOfClosedMOs() ; j++ )
			s << mo.getClosedMO(i, j) << " ";
		s << endl;
	}			
*/
}

/*
istream& operator>>(istream & s, extMOsSel & base)
{
	s >> ((extMOsBase &) base);
	s >> ((Container<extEntry> &) base);
	return s;	
}
*/

void	extMOsSel::threshCut(EnergyType Threshold, INT divideByCSFs)
{
Pix	p = NULL;

	for ( ContainerIterator iter = first() ; !isLast(iter) ;  )
		if ( extEntry	*ext = operator [] (iter) )
		{
		const EnergyType	*dE = ext->getEnergyP();
		INT	n =  ext->getNumberOfRoots();
		INT	CSFs =  ext->getNumberOfCSFs();
		INT	flag = 1;
		
			if ( ext->getSelectionFlags() & 2 )
				flag = 0;
			else
				for ( INT i=0 ; i<n ; i++ )
				{
	//				cout << dE[i] << " " << CSFs << " " <<
	//				 (fabs(dE[i] / (divideByCSFs ? CSFs : 1))) << endl;
					if ( fabs(dE[i] / (divideByCSFs ? CSFs : 1))>=Threshold )
					{
						flag = 0;
						break;
					}
				}

			if ( flag )
			{
				if ( p )
				{
					next(iter);
					del_after(p);
				}
				else
				{
					del_front();
					iter = first();
					continue;
				}
			}
			else
			{
				p = iter.pix;
				next(iter);
			}
		}
}



void	extMOsSel::add(
	const Configuration<MOType> & conf,
	INT	roots,
	const EnergyType *energy,
	INT	CSFs,
	const PTCIType *c,
	INT	SelectionFlags
	)
{
extEntry	*ext = new extEntry(this, roots, energy, CSFs, c, conf, SelectionFlags);
	append(ext);
}


void	extMOsSel::add(
		const Configuration<MOType> &conf, const RootEnergies &r)
{
extEntry	*ext = new extEntry(this, r, conf);
	append(ext);
}


void	extMOsSel::test()
{
}
