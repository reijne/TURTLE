//***********************************************************************
//
//	Name:			RootEnergies.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.04.1997
//
//
//
//
//***********************************************************************

#ifndef __RootEnergies_h
#define __RootEnergies_h

#include "../../../config.h"

#include "EnergyType.h"
#include "PTCIType.h"

#include <iostream>
using std::istream;
using std::ostream;

class RootEnergies {
public:
	RootEnergies();
	RootEnergies(istream &s);
	RootEnergies(INT nRoots, const EnergyType *e, INT nCSFs, const PTCIType *ci);
	RootEnergies(INT nRoots, const EnergyType e);
	RootEnergies(EnergyType e);
	~RootEnergies();

	RootEnergies(const RootEnergies &);
	RootEnergies & operator = (const RootEnergies &);
	
	//----------------------------------------------------------------------

	INT	getNumberOfRoots() const;
	INT	getNumberOfCSFs() const;
		
	const RootEnergies &	getEnergy() const;

	EnergyType	getEnergy(INT rootNo) const;
	const EnergyType	*getEnergyP() const;
	void	setEnergy(INT rootNo, EnergyType e);
	
	PTCIType	getCICoef(INT rootNo, INT CSFNo) const;
	const PTCIType	*getCICoefP() const;
	void	setCICoef(INT rootNo, INT CSFNo, PTCIType e);
	
	void	setStorePTEnergy(INT);
	INT	getStorePTEnergy() const;

	void	setStorePTCoef(INT);
	INT	getStorePTCoef() const;

	//---------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const RootEnergies &);
	
	//---------------------------------------------------------------------

protected:
INT	nRoots;						// number of roots
INT	nCSFs;						// number CSFs 
EnergyType	*e;					// pointer to array of energies
PTCIType		*c;				// pointer to array of estimated ci-vectors

static INT	storePTEnergy;		// flag if PT energies are stored in tree
static INT	storePTCoef;		// flag if PT coefficients are stored in tree
};


#include <string>


inline
RootEnergies::RootEnergies()
{
	nRoots = 0;
	e = NULL;
	nCSFs = 0;
	c = NULL;	
}

inline
RootEnergies::RootEnergies(
	INT _nRoots,
	const EnergyType *_e,
	INT _nCSFs,
	const PTCIType *_c)
{
	nRoots = _nRoots;
	e = NULL;
	c = NULL;
	if ( _e )
	{
		e = new EnergyType[nRoots];
		std::memcpy(e, _e, nRoots*sizeof(EnergyType));
	}
	
	nCSFs = _nCSFs;

	if ( storePTCoef && _c)
	{
		c = new PTCIType[nRoots*nCSFs];
		std::memcpy(c, _c, nRoots*nCSFs*sizeof(PTCIType));
	}
}


inline
RootEnergies::RootEnergies(INT _nRoots, const EnergyType _e)
{
	nRoots = _nRoots;
	e = NULL;
	if ( _e )
	{
		e = new EnergyType[nRoots];
		for ( INT i=0 ; i<nRoots ; i++ )
			e[i] = _e;
	}
	nCSFs = 0;
	c = NULL;	
}

inline
RootEnergies::RootEnergies(EnergyType _e)
{
	nRoots = 1;
	e = new EnergyType[nRoots];
	*e = _e;
	nCSFs = 0;
	c = NULL;	
}


inline
RootEnergies::~RootEnergies()
{
	if ( e )
		delete[] e;
	if ( c )
		delete[] c;
}


inline
INT	RootEnergies::getNumberOfRoots() const
{	return nRoots;	}

inline
INT	RootEnergies::getNumberOfCSFs() const
{	return nCSFs;	}


inline
const RootEnergies &	RootEnergies::getEnergy() const
{	return *this;	}

inline
EnergyType	RootEnergies::getEnergy(INT i) const
{	return e[i];	}

inline
const EnergyType *	RootEnergies::getEnergyP() const
{	return e;	}

inline
void	RootEnergies::setEnergy(INT i, EnergyType _e)
{	e[i] = _e;	}


inline
PTCIType	RootEnergies::getCICoef(INT rootNo, INT CSFNo) const
{	return	c[rootNo*nCSFs + CSFNo];	}

inline
const PTCIType	*RootEnergies::getCICoefP() const
{	return	c;	}

inline
void	RootEnergies::setCICoef(INT rootNo, INT CSFNo, PTCIType _c)
{	c[rootNo*nCSFs + CSFNo] = _c;	}
	


	
inline
void	RootEnergies::setStorePTEnergy(INT i)
{	storePTEnergy = i;	}

inline
INT	RootEnergies::getStorePTEnergy() const
{	return storePTEnergy;	}

inline
void	RootEnergies::setStorePTCoef(INT i)
{	storePTCoef = i;	}

inline
INT	RootEnergies::getStorePTCoef() const
{	return storePTCoef;	}


inline
RootEnergies::RootEnergies(const RootEnergies & re)
{
	e = NULL;
	nRoots = re.nRoots;
	if ( re.e )
	{
		e = new EnergyType[nRoots];
		std::memcpy(e, re.e, nRoots*sizeof(EnergyType));
	}
	c = NULL;
	nCSFs = re.nCSFs;
	if ( re.c )
	{
		c = new PTCIType[nRoots*nCSFs];
		std::memcpy(c, re.c, nRoots*nCSFs*sizeof(PTCIType));
	}
}

inline
RootEnergies &	RootEnergies::operator = (const RootEnergies & re)
{
	if ( e )
		delete[] e;
	e = NULL;
	nRoots = re.nRoots;
	if ( re.e )
	{
		e = new EnergyType[nRoots];
		std::memcpy(e, re.e, nRoots*sizeof(EnergyType));
	}

	if ( c )
		delete[] c;
	c = NULL;
	nCSFs = re.nCSFs;
	if ( re.c )
	{
		c = new PTCIType[nRoots*nCSFs];
		std::memcpy(c, re.c, nRoots*nCSFs*sizeof(PTCIType));
	}
	return *this;
}





#endif
