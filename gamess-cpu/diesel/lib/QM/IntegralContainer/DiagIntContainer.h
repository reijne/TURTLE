#include "../../../config.h"
//***********************************************************************
//
//	Name:			DiagIntContainer.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __DiagIntContainer_h
#define __DiagIntContainer_h


#include "../IO/Fortran/Fort31File.h"
#include "IntegralType.h"
#include "../MRTree/EnergyType.h"

#include "../Symmetry/IrRep.h"

#include "../MO/MOIrReps.h"
#include "../Configuration/Configuration.h"
#include "../Configuration/ConfigurationSet.h"

#include "../../Container/SortedList.h"


template <class MatrixType, class VectorType> class	HMatElements;
class	MOEquivalence;

class DiagIntContainer {
public:
	DiagIntContainer(Fort31File f31, INT Multiplicity);
	~DiagIntContainer();

typedef double	MatrixType;
typedef double	VectorType;

	void	setOneInt(MOType, OneElectronIntegralType);
	void	setTwoIntIIII(MOType, TwoElectronIntegralType);
	void	setTwoIntIIJJ(MOType, MOType, TwoElectronIntegralType);
	void	setTwoIntIJJI(MOType, MOType, TwoElectronIntegralType);

	OneElectronIntegralType	getOneInt(MOType) const;
	TwoElectronIntegralType	getTwoIntIIII(MOType) const;
	TwoElectronIntegralType	getTwoIntIIJJ(MOType, MOType) const;
	TwoElectronIntegralType	getTwoIntIJJI(MOType, MOType) const;

	INT	set(MOType *ind, TwoElectronIntegralType a);



struct	TConf
	{
		TConf() : e(0) {} 
		TConf(EnergyType e, const Configuration<MOType> & conf) :
			e(e), conf(conf) {}
			
		EnergyType				e;
		Configuration<MOType>	conf;
		
		bool	operator < (const TConf &c) const
		{	return e<c.e;	}
	};

	Configuration<MOType>	calcGroundStateClosed(INT nElectrons);
	ConfigurationSet	calcRefConfs(
		Configuration<MOType> groundState, INT nConfs);

	ConfigurationSet	calcRefConfs(
		INT nElectrons,
		IrRep irrep,
		Configuration<MOType> &groundState,
		INT nConfs,
		INT killDegen);

	SortedList<TConf>	*	calcRefConfs(
		INT nElectrons,
		INT nConfs);


	MOEquivalence *	calcMODegeneration(INT nElectrons);

	
	

private:
	void	load1Integrals();
	void	load2Integrals();
	EnergyType	calcDiagClosed(INT *occupation);
	EnergyType	calcDiagClosed(const MOType *mos, INT nMOs) const;
	EnergyType	calcDiagClosed(const MOType *mos, INT nMOs, MOType mo1) const;
	EnergyType	calcDiag(const Configuration<MOType> &conf) const;
	void	dist(MOType nMOsInIrRep, IrRep nSym);
	void	performExcitation(
		INT toAnnihilate, INT toCreate, 
		MOType mo1, MOType mo2, 
		Configuration<MOType> conf);
	void	createActiveSpace(Configuration<MOType> conf);

Fort31File f31;
MOIrReps	*moirreps;
double	core;							//	core enrgy from fort31 file
MOType	maxMO;							//	highest MO
INT	nIrreps;							//	number of irreps
IrRep	irrep;							//	irrep of ground state
INT	*occInd;								//	occupation pattern
INT	*minInd;							//	minimum occupation pattern
INT	nMOs;								//	number of doubly occupied MOs (=nElectrons/2)
EnergyType	minDiag;					//	minimum diagonal energy
OneElectronIntegralType	*pOne;			//	pointer to one electron integrals
TwoElectronIntegralType	*pTwoIIII;		//	pointer to two electron integrals
TwoElectronIntegralType	*pTwoIIJJ;		//	pointer to two electron integrals
TwoElectronIntegralType	*pTwoIJJI;		//	pointer to two electron integrals


INT	n;									//	nth generated configuration
SortedList<TConf>	*lowest;
/*INT	nConfs;								//	size of top list
INT	nInList;							//	actually contained number of 
										//	configurations in top nConfs list
TConf	*confs;							//	top list
*/
INT	minOpenShells;						//	number of minimal open shells
INT	maxOpenShells;						//	number of maximal open shells

HMatElements<MatrixType, VectorType>	**repMatsP5;	//	representation matrices for
										//	diagonal case
										
										
static const INT	activeOccupied = 4;			//	number of highest occupied MOs 
										//	in specific irrep which may be 
										//	annihilated

static const INT	activeVirtual = 4;			//	number of lowest unoccupied MOs 
										//	in specific irrep which may be 
										//	created

INT	*active;							//	flag if specific MO is active

};


inline
void	DiagIntContainer::setOneInt(MOType i, OneElectronIntegralType a)
{	pOne[i-1] = a;	}

inline
void	DiagIntContainer::setTwoIntIIII(MOType i, TwoElectronIntegralType a)
{	pTwoIIII[i-1] = a;	}

inline
void	DiagIntContainer::setTwoIntIIJJ(MOType i, MOType j, TwoElectronIntegralType a)
{
	if ( i>j )
		pTwoIIJJ[(i-1)*maxMO + (j-1)] = a;
	else
		pTwoIIJJ[(j-1)*maxMO + (i-1)] = a;
}

inline
void	DiagIntContainer::setTwoIntIJJI(MOType i, MOType j, TwoElectronIntegralType a)
{
	if ( i>j )
		pTwoIJJI[(i-1)*maxMO + (j-1)] = a;
	else
		pTwoIJJI[(j-1)*maxMO + (i-1)] = a;
}

inline
OneElectronIntegralType	DiagIntContainer::getOneInt(MOType i) const
{
	return pOne[i-1];
}

inline
TwoElectronIntegralType	DiagIntContainer::getTwoIntIIII(MOType i) const
{
	return pTwoIIII[i-1];
}

inline
TwoElectronIntegralType	DiagIntContainer::getTwoIntIIJJ(MOType i, MOType j) const
{
	if ( i>j )
		return pTwoIIJJ[(i-1)*maxMO + (j-1)];
	else
		return pTwoIIJJ[(j-1)*maxMO + (i-1)];
}

inline
TwoElectronIntegralType	DiagIntContainer::getTwoIntIJJI(MOType i, MOType j) const
{
	if ( i>j )
		return pTwoIJJI[(i-1)*maxMO + (j-1)];
	else
		return pTwoIJJI[(j-1)*maxMO + (i-1)];
}


inline
INT	DiagIntContainer::set(MOType *ind, TwoElectronIntegralType a)
{
	if ( ind[0]==ind[1] && ind[0]==ind[2] && ind[0]==ind[3] )
	{
		setTwoIntIIII(ind[0], a);
		return 1;
	}
	if ( ind[0]==ind[1] && ind[2]==ind[3] )
	{
		setTwoIntIIJJ(ind[0], ind[2], a);
		return 1;
	}
	if ( ind[0]==ind[3] && ind[1]==ind[2] || ind[0]==ind[2] && ind[1]==ind[3] )
	{
		setTwoIntIJJI(ind[0], ind[1], a);
		return 1;
	}
	return 0;
}

#endif
