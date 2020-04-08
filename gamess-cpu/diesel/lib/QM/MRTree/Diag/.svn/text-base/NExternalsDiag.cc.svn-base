//***********************************************************************
//
//	Name:			NExternalsDiag.cc
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

#include "../Set/NExternalsSet.h"
#include "../Container/ContainerIterator.h"

#include "../../IO/Array.h"
#include "../../Configuration/Excitation.h"
#include "../../../Math/etc/BinomialCoefficient.h"
#include "../../Configuration/TableCase.h"

#include "../../Cache/VarSizeReadOnlyCache.h"
#include "../../IO/Fortran/Fort31File.h"


#include <stdlib.h>
#include <fstream>
#include "../../../Container/String.h"

#include "../../IO/Verbosity.h"



NExternalsDiag::NExternalsDiag(const NExternalsDiag &t) :
	NExternalsBase<InternalConfsDiag>(t),
        Tree<void>(NULL),
        IndexedContainer<InternalConfsDiag>(t.getNumberOfElements())

{
	for ( ContainerIterator iter = t.first() , iter1 = first() ; 
		!t.isLast(iter) ; 
		t.next(iter), next(iter1) )
	{
		operator [] (iter1) = new InternalConfsDiag(*t[iter]);
		operator [] (iter1)->setParent(this);
	}
	initAll();
}



NExternalsDiag::NExternalsDiag(istream &s, const MOIrReps *moIrReps) :
	NExternalsBase<InternalConfsDiag>(s, moIrReps),
        Tree<void>((void *) NULL),
        IndexedContainer<InternalConfsDiag>(s)
{
	References = 0;
	RefSAFs = 0;
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
	{
	InternalConfsDiag	*p = operator [] (iter);
		p->setParent(this);
		for ( ContainerIterator iter2 = p->first() ; !p->isLast(iter2) ; p->next(iter2) )
			if ( p->operator [] (iter2)->isReference() )
			{
				References++;
				addRefSAFs(getNumberOfSpinAdaptedFunctions(
					p->operator [] (iter2)->getNumberOfOpenShells()));
			}
	}


	initAll();

	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
		cout << "intern-" << i << ": " 
			<< operator [](i)->getNumberOfElements() << endl;

	cout << endl;
	cout << "number of reference configurations          : " << 
		getNumberOfReferences() << endl;	

	cout << "number of CSFs from reference configurations: " << 
		getNumberOfRefConfSpinAdaptedFunctions() << endl;	

	cout << "total number of configurations              : " << 
		getNumberOfLeaves() << endl;	

	cout << "total number of CSFs                        : " << 
		getNumberOfTotalSpinAdaptedFunctions() << endl;	

	cout << endl;
}


/*
NExternalsDiag::NExternalsDiag(MRConfInput &mrconf) :
	NExternalsBase(
		mrconf.getMRMOs(), 
		mrconf.getNumberOfElectrons(), 
		mrconf.getMultiplicity(),
		mrconf.getExcitationLevel()), 
	IndexedContainer<InternalConfsDiag>(mrconf.getExcitationLevel()+1),
	Tree<void>((void *) NULL)
{
	totalSymmetry = mrconf.getRefConf(0).calcIrRep(*mrconf.getMRMOs());
	
InternalConfsDiag	*intern = new InternalConfsDiag (this);
	(*this)[0] = intern;
	for ( INT i=1 ; i<mrconf.getExcitationLevel()+1 ; i++ )
		(*this)[i] = NULL;
		
	RefSAFs = 0;
	for ( INT i=0 ; i<mrconf.getNumberOfRefConfs() ; i++ )
	{
	TupelStructureDiag *tupel = new TupelStructureDiag(
			1, totalSymmetry,
			mrconf.getRefConf(i),
			intern,
			1);
		intern->add(tupel);
	 extMOsDiag	*extMOs = new extMOsDiag(1, 1, 0, 0, tupel);
		(*tupel)[0] = extMOs;
		extMOs->setEnergy(999);
	}
}
*/

NExternalsDiag::NExternalsDiag(const NExternalsSet &set, INT referenceFlag) :
        NExternalsBase<InternalConfsDiag>(set.getMRMOs(), set.getNumberOfElectrons(), set.getMultiplicity()),
        Tree<void>((void *) NULL),
        IndexedContainer<InternalConfsDiag>(set.length())
{
ContainerIterator j=set.first();
	for ( INT i=0 ; 
		i<getNumberOfElements() ; i++ , set.next(j) )
	{
//		cout << "NExt " << i << endl;
		p[i] = new InternalConfsDiag(*set[j], this, referenceFlag);
//		cout << "NExt " << i << endl;
	}

	if ( referenceFlag )
		References = getNumberOfLeaves();
	else
		References = 0;
		
	NumberOfRoots = set.getNumberOfRoots();
	RefSAFs = SAFs = init(0, mrmos->getMOSymmetryP());
}


/*
NExternalsDiag::NExternalsDiag(const ConfigurationSet &set) :
	NExternalsBase(set.getMRMOs(), set.getNumberOfElectrons(), set.getMultiplicity()),
	IndexedContainer<InternalConfsDiag>(set.length()),
	Tree<void>((void *) NULL)
{
ContainerIterator j=set.first();
	for ( INT i=0 ; 
		i<getNumberOfElements() ; i++ , set.next(j) )
	{
//		cout << "NExt " << i << endl;
		p[i] = new InternalConfsDiag(*set[j], this, referenceFlag);
//		cout << "NExt " << i << endl;
	}

	if ( referenceFlag )
		References = getNumberOfLeaves();
	else
		References = 0;
		
	NumberOfRoots = set.getNumberOfRoots();
	RefSAFs = SAFs = init(0, mrmos->getMOSymmetryP());
}

*/


//============================================================================


void	NExternalsDiag::initAll()
{
	NExternalsBase<InternalConfsDiag>::initAll();
	SAFs = init(0, mrmos->getMOSymmetryP());
	
}

NExternalsDiag NExternalsDiag::projectOnReferences() const
{
NExternalsDiag	projected(
	mrmos, electrons, multiplicity, 0);

const InternalConfsDiag	*p = (*this)[first()];
	projected[0] = new InternalConfsDiag(
//		0, p->getNumberOfElements(), (NExternalsDiag *) this);
		0, getNumberOfReferences(), (NExternalsDiag *) this);
INT i = 0;


	for ( ContainerIterator j = p->first() ; !p->isLast(j) ; p->next(j) )
	{
		if ( (*p)[j]->isReference() )
			(*projected[0])[i++] = new TupelStructureDiag(*(*p)[j]);
	}

	return NExternalsDiag(projected);
}


void	NExternalsDiag::writeToStream(ostream & s) const
{
	NExternalsBase<InternalConfsDiag>::writeToStream(s);
	Container<InternalConfsDiag>::writeToStream(s);
}


ostream& operator<<(ostream & s, const NExternalsDiag &base)
{
	s << ((NExternalsBase<InternalConfsDiag> &) base);
	return s;	
}

/*
istream& operator>>(istream & s, NExternalsDiag &base)
{
	s >> ((NExternalsBase &) base);
	s >> ((Container<InternalConfsDiag> &) base);
	return s;	
}
*/






