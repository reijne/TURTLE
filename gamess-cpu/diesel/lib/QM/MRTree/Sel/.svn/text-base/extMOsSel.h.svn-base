//***********************************************************************
//
//	Name:			extMOsSel.h
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

#ifndef __extMOsSel_H
#define __extMOsSel_H

#include "../../../../config.h"

#include "../Base/extMOsBase.h"

#include "../../MO/MOType.h"
#include "../../Configuration/Configuration.h"
#include "../../Configuration/ConfigurationSAFNr.h"

#include "TupelStructureSel.h"
#include "../MRTreeIterator.h"

#include "../Base/extEntry.h"
#include "../Container/ListContainerTree.h"

#include <iostream>

#include "../EnergyType.h"


class	TupelStructureSel;
class	MOEquivalence;


class extMOsSel  : 
	public ListContainerTree<TupelStructureBase, extEntry>,
	public extMOsBase {
public:
	extMOsSel(TupelStructureBase *parent);
	extMOsSel(istream &s):
	Tree<TupelStructureBase>(NULL),
	ListContainer<extEntry>(s),
        extMOsBase(s) {}
	extMOsSel(INT Syms, INT open, 
		INT closed, TupelStructureBase *parent);
/*		 :
		extMOsBase(number, Syms, open, closed),
		Tree<TupelStructureBase>(parent) {}
*/

	~extMOsSel();
	
	extMOsSel(const extMOsSel &);

	//----------------------------------------------------------------

//	MOType	* operator [] (INT n);

	ConfigurationSAFNr<MOType>	getConfigurationSAFNr(INT n) const;
	
	ConfigurationSAFNr<MOType>	
		getConfigurationSAFNr(const ContainerIterator &) const;

	ConfigurationSAFNr<MOType>	
		getConfigurationSAFNr(Pix) const;

	Pix	seek(Configuration<MOType> conf) const;
	
	
//	MOType & operator [] (ContainerIterator) const;


	TupelStructureSel * getParent() const;

	INT calcNumberOfLeaves() const;
	INT getNumberOfLeaves() const;
	
	//----------------------------------------------------------------

	void	setEnergy(const RootEnergies & e);
	void	setEnergy(const ContainerIterator &n, const RootEnergies & e);
	const RootEnergies &	getEnergy(
		const ContainerIterator &n = (INT) 0) const;

	INT	getSelectionFlags(const ContainerIterator &n = (INT) 0) const;


	void	setOpenMO(const ContainerIterator &n, INT i, MOType mo);
	void	setClosedMO(const ContainerIterator &n, INT i, MOType mo);
	
	MOType	getOpenMO(const ContainerIterator &n, INT i) const;
	MOType	getClosedMO(const ContainerIterator &n, INT i) const;

//--------------------------------------------------------------------------	
	
	void	add(
		const Configuration<MOType> & conf,
		INT	roots,
		const EnergyType *energy,
		INT	CSFs,
		const PTCIType *c,
		INT	SelectionFlags);

	void	add(
		const Configuration<MOType> &, const RootEnergies &);

//--------------------------------------------------------------------------	

	void	cutTreeLevel2() {}

	void	threshCut(EnergyType Threshold, INT divideByCSFs);

//----------------------------------------------------------------	

	void	test();

	INT	init(INT, const IrRep *);

//--------------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, extMOsSel &);
	friend istream& operator>>(istream & s, extMOsSel &);

//--------------------------------------------------------------------------

private:
};




inline
TupelStructureSel * extMOsSel::getParent() const
{	return ((TupelStructureSel *) extMOsBase::getParent());	}

		
//----------------------------------------------------------------


inline
INT	extMOsSel::getNumberOfLeaves() const
{	return length();  }

inline
INT	extMOsSel::calcNumberOfLeaves() const
{	return length();  }

//----------------------------------------------------------------

inline
ConfigurationSAFNr<MOType>	extMOsSel::
	getConfigurationSAFNr(const ContainerIterator & iter) const
{	return getConfigurationSAFNr(iter.pix);	}

//----------------------------------------------------------------


inline
void	extMOsSel::setEnergy(const ContainerIterator &n, const RootEnergies & e)
{	operator [] (n)->setEnergy(e);	}


inline
void	extMOsSel::setEnergy(const RootEnergies & e)
{	
	if ( first().pix )
		operator [] (first())->setEnergy(e);
	else
	{
		extEntry	*ext = new extEntry(this, e, 0, NULL);
		append(ext);
	}
}


inline
const RootEnergies &	extMOsSel::getEnergy(const ContainerIterator &n) const
{	return ((extMOsSel *) this)->operator [] (n)->getEnergy();	}

inline
INT	extMOsSel::getSelectionFlags(const ContainerIterator &n) const
{	return ((extMOsSel *) this)->operator [] (n)->getSelectionFlags();	}

//----------------------------------------------------------------

inline
void	extMOsSel::setOpenMO(const ContainerIterator &n, INT i, MOType _mo)
{	operator [] (n)->setMO(i, _mo);	}

inline
void	extMOsSel::setClosedMO(const ContainerIterator &n, INT i, MOType _mo)
{	operator [] (n)->setMO(i+open, _mo);	}

inline
MOType	extMOsSel::getOpenMO(const ContainerIterator &n, INT i) const
{	return ((extMOsSel *) this)->operator [] (n)->getMO(i);	}

inline
MOType	extMOsSel::getClosedMO(const ContainerIterator &n, INT i) const
{	return ((extMOsSel *) this)->operator [] (n)->getMO(i+open);	}


#endif
