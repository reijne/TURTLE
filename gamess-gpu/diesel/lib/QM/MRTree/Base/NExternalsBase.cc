//***********************************************************************
//
//	Name:			NExternalsBase.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.03.1997
//
//
//
//
//
//***********************************************************************

#include <iostream>

#include "NExternalsBase.h"
#include "InternalConfsBase.h"
#include "TupelStructureBase.h"
#include "extMOsBase.h"

#include "../Diag/TupelStructureDiag.h"
#include "../Sel/TupelStructureSel.h"
#include "../Set/TupelStructureSet.h"

#include "../../Configuration/Excitation.h"
#include "../../../Math/etc/BinomialCoefficient.h"
#include "../../Configuration/TableCase.h"

#include "../../Cache/VarSizeReadOnlyCache.h"
#include "../../IO/Fortran/Fort31File.h"

#include <stdlib.h>
#include <fstream>
#include "../../../Container/String.h"

using namespace std;

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::NExternalsBase()
{
	SAFs = 0;
	RefSAFs = 0;
	NumberOfRoots = 0;
	rootNumbers = NULL;
	totalSymmetry = 0;
	mrmos = NULL;
	electrons = 0;
	multiplicity = 0;
	eigFuncs = NULL;
        References = 0;
}

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::
	NExternalsBase(
	const MOIrReps *_moirreps,
	INT NumberOfElectrons, INT Multiplicity)
	
{
	SAFs = 0;
	RefSAFs = 0;
	if ( _moirreps )
	{
		mrmos = new MRMOs;
		*mrmos = *_moirreps;
	}
	totalSymmetry = 0;
	electrons = NumberOfElectrons;
	if ( electrons>MAXELECTRONS )
	{	cout << "error: number of electrons too large" << endl;
		cout << "Please increase parameter MAXELECTRONS in file ConfigurationGlobals.h" << endl;
		cout << "and recompile." << endl;
		exit(1);
	}
	multiplicity = Multiplicity;
//	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
//		(*this)[i] = NULL;
	eigFuncs = new SpinEigenFunctionDegeneration(multiplicity, MAXOPENSHELLS);
	NumberOfRoots = 0;
	rootNumbers = NULL;
        References = 0;
}

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::
	NExternalsBase(
	const MRMOs *_mrmos,
	INT NumberOfElectrons, INT Multiplicity)
	
{
	SAFs = 0;
	RefSAFs = 0;
	if ( _mrmos )
	{
		mrmos = new MRMOs;
		*mrmos = *_mrmos;
	}
	totalSymmetry = 0;
	electrons = NumberOfElectrons;
	if ( electrons>MAXELECTRONS )
	{	cout << "error: number of electrons too large" << endl;
		cout << "Please increase parameter MAXELECTRONS in file ConfigurationGlobals.h" << endl;
		cout << "and recompile." << endl;
		exit(1);
	}
	multiplicity = Multiplicity;
//	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
//		(*this)[i] = NULL;
	eigFuncs = new SpinEigenFunctionDegeneration(multiplicity, MAXOPENSHELLS);
	NumberOfRoots = 0;
	rootNumbers = NULL;
        References = 0;
}

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::
	~NExternalsBase()
{
	if ( eigFuncs )
		delete eigFuncs;
//	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
//		delete (*this)[i];
	if ( mrmos )
		delete mrmos;
	if ( rootNumbers )
		delete[] rootNumbers;
}

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::NExternalsBase(istream &s, const MOIrReps *_moIrReps)
{
//	cout << "NExternalsBase(istream &s)" << endl;
	if ( _moIrReps )
		mrmos = new MRMOs(*_moIrReps);
	else
		mrmos = NULL;
	s >> SAFs;
	s >> RefSAFs;
	s >> NumberOfRoots;
	rootNumbers = new INT[NumberOfRoots];
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		s >> rootNumbers[i];
	s >> totalSymmetry;
	s >> electrons;
	s >> multiplicity;
	s >> References;
	eigFuncs = new SpinEigenFunctionDegeneration(multiplicity, MAXOPENSHELLS);
//	cout << "Test:---------------------" << endl;
//	cout << *this << endl;
}




template <class ContainedObjectType>
NExternalsBase<ContainedObjectType>::NExternalsBase(const NExternalsBase<ContainedObjectType> & nExt)
{
	SAFs = nExt.getNumberOfTotalSpinAdaptedFunctions();
	RefSAFs = nExt.getNumberOfRefConfSpinAdaptedFunctions();
	NumberOfRoots = nExt.getNumberOfRoots();
	rootNumbers = new INT[NumberOfRoots];
	memcpy(rootNumbers, nExt.getRootNumbersP(), NumberOfRoots*sizeof(INT));
	totalSymmetry = nExt.getTotalSymmetry();
	electrons = nExt.getNumberOfElectrons();
	multiplicity = nExt.getMultiplicity();
	References = nExt.getNumberOfReferences();
	
	if ( nExt.mrmos )
	{
		mrmos = new MRMOs;
		*mrmos = *nExt.mrmos;
	}
	else
		mrmos = NULL;
	eigFuncs = new SpinEigenFunctionDegeneration(*nExt.eigFuncs);
	internalIntersection = nExt.internalIntersection;
}

template <class ContainedObjectType>
NExternalsBase<ContainedObjectType> &	NExternalsBase<ContainedObjectType>::operator = (const NExternalsBase<ContainedObjectType> & nExt)
{
	if ( eigFuncs )
		delete eigFuncs;
	if ( mrmos )
		delete mrmos;
	if ( rootNumbers )
		delete rootNumbers;

	SAFs = nExt.getNumberOfTotalSpinAdaptedFunctions();
	RefSAFs = nExt.getNumberOfRefConfSpinAdaptedFunctions();
	NumberOfRoots = nExt.getNumberOfRoots();
	rootNumbers = new INT[NumberOfRoots];
	memcpy(rootNumbers, nExt.getRootNumbersP(), NumberOfRoots*sizeof(INT));
	totalSymmetry = nExt.getTotalSymmetry();
	electrons = nExt.getNumberOfElectrons();
	multiplicity = nExt.getMultiplicity();
	References = nExt.getNumberOfReferences();
	
	mrmos = NULL;
	eigFuncs = NULL;
	internalIntersection = nExt.internalIntersection;
	return *this;
}


template <class ContainedObjectType>
void	NExternalsBase<ContainedObjectType>::calcInternalIntersection()
{
INT	f = 1;
Configuration<MOType>	conf;
	for ( ContainerIterator i = this->first() ; !this->isLast(i) ; this->next(i) )
	{
	ContainedObjectType	*p = (*this)[i];
		for ( ContainerIterator j = p->first() ; !p->isLast(j) ; p->next(j) )
		{
			if ( f )
				conf = (Configuration<MOType> &) *(*p)[j];
			else
				conf &= (Configuration<MOType> &) *(*p)[j];
			f = 0;
		}
	}
	cout << endl;
	cout << "internal intersection: " << conf << endl;
	cout << endl;
}

template <class ContainedObjectType>
ConfigurationSet	NExternalsBase<ContainedObjectType>::getRefConfSet() const
{
ConfigurationSet	set;
	for ( ContainerIterator i = this->first() ; !this->isLast(i) ; this->next(i) )
	{
	const ContainedObjectType	*p = (*this)[i];
		for ( ContainerIterator j = p->first() ; !p->isLast(j) ; p->next(j) )
			if ( (*p)[j]->isReference() )
			{
			Configuration<MOType>	conf = (Configuration<MOType> &) *(*p)[j];
				set.add(conf);
			}
	}
	return set;
}



template <class ContainedObjectType>
void	NExternalsBase<ContainedObjectType>::initAll()
{
//	for ( INT ii=0 ; ii<getNumberOfLeaves() ; ii++ )
//		cout << ii << ": " << (*this)(ii) << endl;


	
/*	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
		for ( INT j=0 ; j<(*this)[i]->getNumberOfElements() ; j++ )
			for ( INT k=0 ; k<(*(*this)[i])[j]->getNumberOfElements() ; k++ )
			{//	printf("%d %d %d\n", i, j, k);
				(*(*(*this)[i])[j])[k]->init();
			}
*/

/*	DirPath += "/fort.31";
Fort31File	f31;
	f31.name = DirPath.chars();
	f31.format = format;
*/
INT	first = 1;
Configuration<MOType>	inactive;


	for ( INT i=0 ; i<this->getNumberOfElements() ; i++ )
		for ( INT j=0 ; j<(*this)[i]->getNumberOfElements() ; j++ )
		{
			if ( (*(*this)[i])[j]->isReference() )
			{	for ( INT k=0 ; k<(*(*this)[i])[j]->getNumberOfOpenShells() ; k++ )
					mrmos->setInternal((*(*this)[i])[j]->getOpenShell(k));
				for ( INT k=0 ; k<(*(*this)[i])[j]->getNumberOfClosedShells() ; k++ )
					mrmos->setInternal((*(*this)[i])[j]->getClosedShell(k));
				if ( first )
				{
					inactive = *(*(*this)[i])[j];
					first = 0;
				}
				else
					inactive &= *(*(*this)[i])[j];
			}
		}
	
	mrmos->setNumberOfInactiveMOs(
		0*inactive.getNumberOfOpenShells() + inactive.getNumberOfClosedShells());

/*	for ( INT i=0 ; i<inactive.getNumberOfOpenShells() ; i++ )
		mrmos->setInactiveMO(i, inactive.getOpenShell(i));
*/
	for ( INT i=0 ; i<inactive.getNumberOfClosedShells() ; i++ )
		mrmos->setInactiveMO(0*inactive.getNumberOfOpenShells() + i,
			 inactive.getClosedShell(i));

	mrmos->initIntExt();
}





template <class ContainedObjectType>
void	NExternalsBase<ContainedObjectType>::writeToStream(ostream & s) const
{
	s << SAFs << endl;
	s << RefSAFs << endl;
	s << NumberOfRoots << endl;
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		s << rootNumbers[i] << endl;
	s << totalSymmetry << endl;
	s << electrons << endl;
	s << multiplicity << endl;
	s << References << endl;
}

template <class ContainedObjectType>
ostream& operator<<(ostream & s, const NExternalsBase<ContainedObjectType> &base)
{
	s << base.SAFs << endl;
	s << base.RefSAFs << endl;
	s << base.NumberOfRoots << endl;
	for ( INT i=0 ; i<base.NumberOfRoots ; i++ )
		s << base.rootNumbers[i] << endl;
	s << base.totalSymmetry << endl;
	s << base.electrons << endl;
	s << base.multiplicity << endl;
	s << base.References << endl;
	return s;	
}

template <class ContainedObjectType>
istream& operator>>(istream & s, NExternalsBase<ContainedObjectType> &base)
{
	s >> base.SAFs;
	s >> base.RefSAFs;
	s >> base.totalSymmetry;
	s >> base.electrons;
	s >> base.multiplicity;
	s >> base.References;
	return s;	
}


#include "../Diag/InternalConfsDiag.h"
#include "../Sel/InternalConfsSel.h"
#include "../Set/InternalConfsSet.h"


template class NExternalsBase<InternalConfsDiag>;
template class NExternalsBase<InternalConfsSel>;
template class NExternalsBase<InternalConfsSet>;

template ostream& operator << (ostream&, const NExternalsBase<InternalConfsDiag> &);
template istream& operator>>(istream & s, NExternalsBase<InternalConfsDiag> &);
template ostream& operator << (ostream&, const NExternalsBase<InternalConfsSel> &);
template istream& operator>>(istream & s, NExternalsBase<InternalConfsSel> &);
template ostream& operator << (ostream&, const NExternalsBase<InternalConfsSet> &);
template istream& operator>>(istream & s, NExternalsBase<InternalConfsSet> &);

