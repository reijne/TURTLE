//***********************************************************************
//
//	Name:	ConfigurationSAFNr.h
//
//	Description:	derived from Configuration
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			06.11.1996
//
//
//
//
//
//***********************************************************************


#ifndef __CONFIGURATIONSAFNR_H
#define __CONFIGURATIONSAFNR_H

#include "../../../config.h"


#include "Configuration.h"
#include "../MRTree/RootEnergies.h"

template <class TMOType> class ConfigurationSAFNr;
template <class TMOType> ostream& operator<< (ostream & s, const ConfigurationSAFNr<TMOType> &);
template <class TMOType> istream& operator>> (istream & s, ConfigurationSAFNr<TMOType> &);


template <class TMOType>
class ConfigurationSAFNr : public Configuration<TMOType> {
public:
	ConfigurationSAFNr(
		Configuration<TMOType> & conf,
		Configuration<TMOType> & internal,
		INT	nExtOpen, INT nExtClosed,
		INT SAFNr, INT SAFInc,
		const RootEnergies &rootEnergies,
		INT isRef,
		INT	SelectionFlags);
	
	ConfigurationSAFNr(
		Configuration<TMOType> & conf,
		Configuration<TMOType> & internal,
		INT	nExtOpen, INT nExtClosed,
		INT SAFNr, INT SAFInc,
		INT roots,
		const EnergyType *Energies,
		INT CSFs,
		const PTCIType *ci,
		INT isRef,
		INT	SelectionFlags);
	
	
	void	setSAFNr(INT i);
	INT		getSAFNr() const;

	void	setSAFInc(INT i);
	INT		getSAFInc() const;
	
	void	setEnergy(RootEnergies e);
	const RootEnergies &	getEnergy() const;
	
	INT	getNExt() const;
	INT	getNExtOpen() const;
	INT	getNExtClosed() const;

	INT	isReference() const;
	INT	getSelectionFlags() const;

	void	setInternal(Configuration<TMOType>);
	Configuration<TMOType>	getInternal() const;
	
//--------------------------------------------------------------------------

	friend ostream& operator<< <TMOType> (ostream & s, const ConfigurationSAFNr<TMOType> &);
	friend istream& operator>> <TMOType> (istream & s, ConfigurationSAFNr<TMOType> &);

//--------------------------------------------------------------------------

	
private:
INT	SAFNr;							// start number of SAF
INT	SAFInc;							// SAFs per configuration
RootEnergies	rootEnergies;		// perturbation energy for certain roots
INT	nExtOpen;						// number of open external MOs
INT nExtClosed;						// number of closed external MOs
Configuration<TMOType>	internal;	// internal part of configuration
INT	isRef;							// flag if reference configuration
INT	SelectionFlags;						// flag if near generation
};

template <class TMOType>	ostream& operator<<(ostream & s, const ConfigurationSAFNr<TMOType> &);
template <class TMOType>	istream& operator>>(istream & s, ConfigurationSAFNr<TMOType> &);






#endif
