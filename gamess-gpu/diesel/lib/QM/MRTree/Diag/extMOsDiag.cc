//***********************************************************************
//
//	Name:			extMOsDiag.cc
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

#include "NExternalsDiag.h"
#include "InternalConfsDiag.h"
#include "TupelStructureDiag.h"
#include "extMOsDiag.h"
#include "../Set/extMOsSet.h"

#include <stdlib.h>
#include <string>
#include <iomanip>

#include "../Sel/extMOsSel.h"

#include "../../MO/MOMapping.h"

using namespace std;

static int	cmp(const void *_p1, const void *_p2)
{
const MOType *p1 =(const MOType *) _p1;
const MOType *p2 =(const MOType *) _p2;
	if ( *p1<*p2 )
		return -1;
	if ( *p1>*p2 )
		return 1;
	return 0;
}



extMOsDiag::extMOsDiag(istream &s) :
	Tree<TupelStructureBase>(NULL)
{
	mo = NULL;
	MOInd = NULL;
	IndTrn = NULL;
	ci = NULL;
	energy = NULL;
	SelectionFlags = NULL;

//	cout << "=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?= A" << endl;
        extMOsSel extSel(s);
//	cout << "=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?= B" << endl;
	constructFromExtMOsSel(extSel);
//	cout << "=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?= C" << endl;
//	cout << extSel << endl;


//	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
//		(*((ContainerTree *) this))[iter]->setParent(this);
}


extMOsDiag::extMOsDiag(const extMOsDiag &e) :
        Tree<TupelStructureBase>(NULL),
	extMOsBase(e)
{
	number = e.number;
	roots = e.roots;
	CSFs = e.CSFs;
	ci = NULL;
	MOInd = NULL;
	IndTrn = NULL;
	SelectionFlags = NULL;
	mo = NULL;
	energy = NULL;
	

	if ( e.mo )
	{
		mo = new MOType * [total];
		for ( INT i=0 ; i<total ; i++ )
		{
			mo[i] = new MOType[number*total];
			memcpy(mo[i], e.mo[i], number*total*sizeof(MOType));
		}
	}
	
	if ( e.MOInd )
	{
		MOInd = new INT * [total];
		for ( INT i=0 ; i<total ; i++ )
		{
			MOInd[i] = new INT[maxMO+1];
			memcpy(MOInd[i], e.MOInd[i], (maxMO+1)*sizeof(INT));
		}
	}
	
	if ( e.IndTrn )
	{
		IndTrn = new INT * [total];
		for ( INT i=0 ; i<total ; i++ )
		{
			IndTrn[i] = new INT [number];
			memcpy(IndTrn[i], e.IndTrn[i], number*sizeof(INT));
		}
	}
	
	if ( e.energy )
	{
		energy = new EnergyType [number];
		memcpy(energy, e.energy, number*sizeof(EnergyType));
	}
	
	if ( e.ci )
	{
		ci = new PTCIType[number*roots*CSFs];
		memcpy(ci, e.ci, number*roots*CSFs*sizeof(PTCIType));
	}
	
	if ( e.SelectionFlags )
	{
		SelectionFlags = new char [number];
		memcpy(SelectionFlags, e.SelectionFlags, number*sizeof(char));
	}
	
}


extMOsDiag::extMOsDiag(INT _number, INT _Syms, INT _open, INT _closed, 
	TupelStructureDiag *parent) :
        Tree<TupelStructureBase>(parent),
	extMOsBase(_Syms, _open, _closed)
{
	number = _number;
	ci = NULL;
	MOInd = NULL;
	IndTrn = NULL;
	SelectionFlags = NULL;
	mo = NULL;
	energy = NULL;
	if ( total )
	{
		mo = new MOType * [total];
		for ( INT i=0 ; i<total ; i++ )
			mo[i] = new MOType[number*total];
	}
	if ( number )
		energy = new EnergyType [number];
}

extMOsDiag::~extMOsDiag()
{
	if ( mo )
	{
		for ( INT i=0 ; i<total ; i++ )
			delete[] mo[i];
		delete[] mo;
	}
	if ( MOInd )
	{	for ( INT i=0 ; i<total ; i++ )
			delete[] MOInd[i];
		delete[] MOInd;
	}

	if ( IndTrn )
	{	for ( INT i=0 ; i<total ; i++ )
			delete[] IndTrn[i];
		delete[] IndTrn;
	}

	
	if ( energy )
		delete[] energy;

	if ( ci )
		delete[] ci;

	if ( SelectionFlags )
		delete[] SelectionFlags;
}

extMOsDiag::extMOsDiag(extMOsSet &set, TupelStructureDiag *parent) :
	Tree<TupelStructureBase>(parent)

{
	mo = NULL;
	MOInd = NULL;
	IndTrn = NULL;
	ci = NULL;
	energy = NULL;
	SelectionFlags = NULL;
	constructFromExtMOsSet(set);
}


extMOsDiag::extMOsDiag(extMOsSel &extSel) :
	Tree<TupelStructureBase>(NULL)
{
	mo = NULL;
	MOInd = NULL;
	IndTrn = NULL;
	ci = NULL;
	energy = NULL;
	SelectionFlags = NULL;
	constructFromExtMOsSel(extSel);
}

void	extMOsDiag::constructFromExtMOsSet(extMOsSet &extSet)
{
	(extMOsBase &) *this = (extMOsBase &) extSet;
	number = extSet.getNumberOfElements();
	if ( !number )
		return;
	roots = extSet.operator [] (extSet.first())->getNumberOfRoots();
	CSFs = extSet.operator [] (extSet.first())->getNumberOfCSFs();


	if ( total )
	{
		mo = new MOType * [total];
		for ( INT i=0 ; i<total ; i++ )
			mo[i] = new MOType[number*total];
	}
	
	energy = NULL;
	ci = NULL;

	if (  extSet.operator [] (extSet.first())->getEnergyP() )
		energy = new EnergyType[number*roots];
	if ( extSet.operator [] (extSet.first())->getCICoefP() )
		ci = new PTCIType[number*roots*CSFs];

	SelectionFlags = new char[number];
        for ( int i=0; i<number; i++ )
            SelectionFlags[i] = 0;
	
INT	iRoot = 0;
INT	iCI = 0;
INT	iMO = 0;
INT	j = 0;
	for ( ContainerIterator iter = extSet.first() ; !extSet.isLast(iter) ; 
		extSet.next(iter) )
	{
	const extEntry	*entry = extSet.operator [] (iter);
		if ( total )
			memcpy(mo[0]+iMO, entry->getMOP(), total*sizeof(MOType));
		iMO += total;

		if ( energy )
			memcpy(energy+iRoot, entry->getEnergyP(), roots*sizeof(EnergyType));
		if ( ci )
			memcpy(ci+iCI, entry->getCICoefP(), roots*CSFs*sizeof(PTCIType));
		SelectionFlags[j++] = entry->getSelectionFlags();
		iRoot += roots;
		iCI += roots*CSFs;
	}
}


void	extMOsDiag::constructFromExtMOsSel(extMOsSel &extSel)
{
//	(extMOsBase &) *this = dynamic_cast<extMOsBase &>(const_cast<extMOsSel&>(extSel));
	(extMOsBase &) *this = (extMOsBase &) extSel;
	
//	cout << ((extMOsBase &) extSel) << endl;
//	cout << ((extMOsBase &) *this) << endl;

	number = extSel.getNumberOfElements();
	if ( !number )
		return;
	roots = extSel.operator [] (extSel.first())->getNumberOfRoots();
	CSFs = extSel.operator [] (extSel.first())->getNumberOfCSFs();

//	cout << "number=" << number << endl;
//	cout << "total=" << total << endl;
//	cout << "roots=" << roots << endl;

	if ( total )
	{
		mo = new MOType * [total];
		for ( INT i=0 ; i<total ; i++ )
			mo[i] = new MOType[number*total];
	}

	energy = NULL;
	ci = NULL;

	if (  extSel.operator [] (extSel.first())->getEnergyP() )
		energy = new EnergyType[number*roots];
	if ( extSel.operator [] (extSel.first())->getCICoefP() )
		ci = new PTCIType[number*roots*CSFs];

	SelectionFlags = new char[number];
	
INT	iRoot = 0;
INT	iCI = 0;
INT	iMO = 0;
INT	j = 0;
	for ( ContainerIterator iter = extSel.first() ; !extSel.isLast(iter) ; 
		extSel.next(iter) )
	{
	const extEntry	*entry = extSel.operator [] (iter);
//		printf("iter=%x\n", iter.i);
//		cout << "iMO=" << iMO << endl;
		if ( total )
			memcpy(mo[0]+iMO, entry->getMOP(), total*sizeof(MOType));
		iMO += total;
//		cout << "iRoot=" << iRoot << endl;
//		cout << "iCI=" << iCI << endl;
		if ( energy )
			memcpy(energy+iRoot, entry->getEnergyP(), roots*sizeof(EnergyType));
		if ( ci )
			memcpy(ci+iCI, entry->getCICoefP(), roots*CSFs*sizeof(PTCIType));
		SelectionFlags[j++] = entry->getSelectionFlags();
		iRoot += roots;
		iCI += roots*CSFs;
	}
}

TupelStructureDiag * extMOsDiag::getParent() const
{	
//	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	Bug in gcc 2.7.2:
//	class pointer arithmetic wrong for virtual bases classes if
//	used in inlined code
//	e.g.:
//	(TupelStructureDiag *) (((char *) extMOsBase::getParent())
//	 - sizeof(DiagTreeBase<InternalConfsDiag, extMOsDiag>)));*/
//	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	return ((TupelStructureDiag *) extMOsBase::getParent());
}
		

ConfigurationSAFNr<MOType>	extMOsDiag::getConfigurationSAFNr(INT n) const
{
Configuration<MOType>	conf = *getParent();

//	cout << "n= " << n << ", this = " << this << endl;
	for ( INT i=0 ; i<getNumberOfOpenMOs() ; i++ )
		conf.create(getOpenMO(n, i));

	for ( INT i=0 ; i<getNumberOfClosedMOs() ; i++ )
	{	conf.create(getClosedMO(n, i));
		conf.create(getClosedMO(n, i));
	}


//	cout << "CSFs=" << CSFs << endl;

	return ConfigurationSAFNr<MOType>(
		conf, *getParent(),
		getNumberOfOpenMOs(), getNumberOfClosedMOs(),
		SAFStart + n*SAFInc, SAFInc,
		roots, (energy ? energy+n*roots : (const double *) NULL),
		CSFs, (ci ? ci + n*roots*CSFs : (const double *) NULL),
		getParent()->isReference(),
		SelectionFlags[n]);
}


static int	cmpfunc1(const void *_p1, const void *_p2)
{
const INT *p1 = (const INT *) _p1;
const INT *p2 = (const INT *) _p2;
	if ( p1[0]>p2[0] )
		return 1;
	if ( p1[0]<p2[0] )
		return -1;
	if ( p1[1]>p2[1] )
		return 1;
	if ( p1[1]<p2[1] )
		return -1;
	return 0;
}

static int	cmpfunc2(const INT *p1, const INT *p2)
{
	if ( p1[1]>p2[1] )
		return 1;
	if ( p1[1]<p2[1] )
		return -1;
	if ( p1[0]>p2[0] )
		return 1;
	if ( p1[0]<p2[0] )
		return -1;
	return 0;
}



INT	extMOsDiag::init(INT n, const IrRep *irrep)
{
//	allocate space for help array to generate index translation table
INT	*moh = new INT [(total+1) * number];
	
//	allocate space for pointers to and index translation table
	if ( total>0 )
	{	IndTrn = new INT * [total];
		for ( INT i=0 ; i<total ; i++ )
			IndTrn[i] = new INT [number];
	}

/* FD gives errors and is not used anywhere
int	(*cmpfunc[2])(const INT *, const INT *);

	cmpfunc[0] = cmpfunc1;
	cmpfunc[1] = cmpfunc2;
*/        
	
	char	*SelectionFlagsTemp = new char[number];
		memcpy(SelectionFlagsTemp, SelectionFlags, number*sizeof(char));

	EnergyType	*energyTemp = NULL;
		if ( energy )
		{
			energyTemp = new EnergyType[number*roots];
			memcpy(energyTemp, energy, number*roots*sizeof(EnergyType));
		}

	PTCIType	*ciTemp = NULL;
		if ( ci )
		{
			ciTemp= new PTCIType[number*roots*CSFs];
			memcpy(ciTemp, ci, number*roots*CSFs*sizeof(PTCIType));
		}
					
//	initialize arrays sorted by n-th MO
	for ( INT i=0 ; i<total ; i++ )
	{
/*		for ( INT j=0 ; j<number ; j++ )
		{	for ( INT k=0 ; k<total ; k++ )
				moh[j*(total+1) + k] = mo[0][j*total + k];
			moh[j*(total+1) + total] = j;
		}
*/
		for ( INT j=0 ; j<number ; j++ )
		{	moh[j*(total+1)] = mo[0][j*total + i];
		INT	kk = 0;
			for ( INT k=0 ; k<total ; k++ )
				if ( k!=i )
					moh[j*(total+1) + ++kk] = mo[0][j*total + k];

			moh[j*(total+1) + total] = j;
		}

		qsort(moh, number, (total+1)*sizeof(INT), cmpfunc1);

		for ( INT j=0 ; j<number ; j++ )
		{	mo[i][j*total + i] = moh[j*(total+1)];
		INT	kk = 0;

			for ( INT k=0 ; k<total ; k++ )
				if ( k!=i )
					mo[i][j*total + k] = moh[j*(total+1) + ++kk];

		INT moInd = moh[j*(total+1) + total];		
			if ( i>0 )
				IndTrn[i][j] = moInd;
			else
			{
				SelectionFlags[j] = SelectionFlagsTemp[moInd];

				if ( energy )
					memcpy(energy + j*roots, energyTemp + moInd*roots,
						roots*sizeof(EnergyType));

				if ( ci )
					memcpy(ci + j*roots*CSFs, ciTemp + moInd*roots*CSFs,
						roots*CSFs*sizeof(EnergyType));

			// inverse mapping
				IndTrn[i][moInd] = j;
			}
				
/*			cout << i << ", " << j << ": ";
			for ( INT ll=0 ; ll<total ; ll++ )
				cout << mo[i][j*total+ll] << " ";
			cout << " | " << IndTrn[i][j] << endl;
*/		}
	}
	if ( ciTemp )
		delete[] ciTemp;
	if ( energyTemp )
		delete[] energyTemp;
	delete[] SelectionFlagsTemp;
	delete[] moh;

//	get maximum used MO number
	maxMO = 0;
	for ( INT i=0 ; i<number ; i++ )
		for ( INT j=0 ; j<total ; j++ )
			if ( maxMO<getOpenMO(i, j) )
				maxMO = getOpenMO(i, j);


//	allocate space for pointers to and access lists
	if ( total )
	{
		MOInd = new INT * [total];
		for ( INT i=0 ; i<total ; i++ )
			MOInd[i] = new INT [maxMO+1];
	}
	
//	initialize access lists
	for ( INT i=0 ; i<total ; i++ )
	{	for ( INT j=0 ; j<=maxMO ; j++ )
			MOInd[i][j] = -1;
			
	INT	lastMO = -1;
		for ( INT j=0 ; j<number ; j++ )
			if ( lastMO<getOpenMO(j, i, i) )
			{//	cout << i << " " << j << " " << getOpenMO(j, i, i) << endl;
				MOInd[i][getOpenMO(j, i, i)] = j;
				lastMO = getOpenMO(j, i, i);
			}
//		cout << "===========" << endl;
	}	



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// initialize inSymN-list
//	printf("%d %x\n", Syms, inSymN);
	delete[] inSymN;
	Syms = 8;
	inSymN = new INT[Syms];
//	printf("%x\n", inSymN);
	memset(inSymN, 0, Syms*sizeof(INT));
	if ( total ) 
		for ( INT i=0 ; i<number ; i++ )
		{
		INT	j = getOpenMO(i, 0);
			if ( j>0 && j<=maxMO ) 
				inSymN[irrep[j]]++;
		}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	// initialize SAFStart
	SAFStart = n;
	SAFInc = getParent()->getParent()->getParent()->
			getNumberOfSpinAdaptedFunctions(
			open + getParent()->getNumberOfOpenShells());
	return n + number*SAFInc;
}



void	extMOsDiag::test()
{
	for ( INT i=0 ; i<number ; i++ )
		cout << MOInd[0][i] << " ";
	cout << endl;
}

/*
static int	cmpfunc2First(const MOType *p1, const MOType *p2)
{
	if ( p1[cmpMode]>p2[cmpMode] )
		return 1;
	if ( p1[cmpMode]<p2[cmpMode] )
		return -1;
	return 0;
}
*/

/*
INT	extMOsDiag::findOpenMO(INT i, MOType MO) const
{
	return MOInd[i][MO];
	cout << "Start!" << i << " " << MO << endl;
	cmpMode = i;
MOType	key[2];
	key[i] = MO;
MOType	*p = (MOType *) bsearch(key, mo[i], number,
		total*sizeof(MOType), cmpfunc2First);
		
	cout << p << endl;
	if ( !p )
		return 1<<30;
		
	cout << " ! " << p[0] << " " << p[1] << endl;
		
	return p-mo[i];
}
*/

void	extMOsDiag::writeToStream(ostream & s) const
{
	s << getNumberOfElements() << endl;

const EnergyType	*pr = energy;
const MOType	*pMO = mo[0];
	
MOType* m = new MOType[total];



	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
	{
		s << roots << endl;
		for ( INT j=0 ; j<roots ; j++ )
			s << setiosflags(ios::scientific) << *pr++ << " ";
		if ( roots )
			s << endl;
			
		s << total << endl;

		for ( INT j=0 ; j<total ; j++ )
			m[j] = moMapping.getReal(*pMO++);
			
		qsort(m, open, sizeof(MOType), cmp);
		qsort(m+open, total-open, sizeof(MOType), cmp);
			
			
			
		for ( INT j=0 ; j<total ; j++ )
			s << m[j] << " ";
		if ( total )
			s << endl;
	}
	extMOsBase::writeToStream(s);
	delete[] m;
}

ostream& operator<<(ostream & s, const extMOsDiag & mo)
{
	for ( INT i=0 ; i<mo.getNumberOfConfs() ; i++ )
	{
		{
		MOType* m = new MOType[mo.getNumberOfOpenMOs()];
			for ( INT j=0 ; j<mo.getNumberOfOpenMOs() ; j++ )
				m[j] = moMapping.getReal(mo.getOpenMO(i, j));

			qsort(m, mo.getNumberOfOpenMOs(), sizeof(MOType), cmp);

			for ( INT j=0 ; j<mo.getNumberOfOpenMOs() ; j++ )
				s << m[j] << " ";
			delete[] m;
		}

		if ( mo.getNumberOfOpenMOs() && mo.getNumberOfClosedMOs() )
			s << "| ";

		{
		MOType* m = new MOType[mo.getNumberOfClosedMOs()];

			for ( INT j=0 ; j<mo.getNumberOfClosedMOs() ; j++ )
				m[j] = moMapping.getReal(mo.getClosedMO(i, j));

			qsort(m, mo.getNumberOfClosedMOs(), sizeof(MOType), cmp);

			for ( INT j=0 ; j<mo.getNumberOfClosedMOs() ; j++ )
				s << m[j] << " ";
			delete[] m;
		}
		s << endl;
	}
	return s;
}

/*
istream& operator>>(istream & s, extMOsDiag & mo)
{
	for ( INT i=0 ; i<mo.getNumberOfConfs() ; i++ )
	{	for ( INT j=0 ; j<mo.getNumberOfOpenMOs() ; j++ )
			s >> mo.getOpenMO(i, j);
		for ( INT j=0 ; j<mo.getNumberOfClosedMOs() ; j++ )
			s >> mo.getClosedMO(i, j);
	}			
}
*/

