//***********************************************************************
//
//	Name:	ConfigurationSAFNr.cc
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


#include "ConfigurationSAFNr.h"
#include <iomanip>

using namespace std;

template <class TMOType>
ConfigurationSAFNr<TMOType>::ConfigurationSAFNr(
	Configuration<TMOType> & conf,
	Configuration<TMOType> & _internal,
	INT	_nExtOpen, INT _nExtClosed,
	INT _SAFNr, INT _SAFInc,
	const RootEnergies &_rootEnergies,
	INT	_isRef,
	INT	_SelectionFlags) :
	Configuration<TMOType>(conf)
{
	SAFNr = _SAFNr;
	SAFInc = _SAFInc;
	rootEnergies = _rootEnergies;
	internal = _internal;
	nExtOpen = _nExtOpen;
	nExtClosed = _nExtClosed;
	isRef = _isRef;
	SelectionFlags = _SelectionFlags;
}


template <class TMOType>
ConfigurationSAFNr<TMOType>::ConfigurationSAFNr(
	Configuration<TMOType> & conf,
	Configuration<TMOType> & _internal,
	INT	_nExtOpen, INT _nExtClosed,
	INT _SAFNr, INT _SAFInc,
	INT roots,
	const EnergyType *Energies,
	INT CSFs,
	const PTCIType *ci,
	INT	_isRef,
	INT	_SelectionFlags) :
	Configuration<TMOType>(conf)
{
	SAFNr = _SAFNr;
	SAFInc = _SAFInc;
	rootEnergies = RootEnergies(roots, Energies, CSFs, ci);
	internal = _internal;
	nExtOpen = _nExtOpen;
	nExtClosed = _nExtClosed;
	isRef = _isRef;
	SelectionFlags = _SelectionFlags;
}

template <class TMOType>
void	ConfigurationSAFNr<TMOType>::setSAFNr(INT i)
{	SAFNr = i;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getSAFNr() const
{	return	SAFNr;	}

template <class TMOType>
void	ConfigurationSAFNr<TMOType>::setSAFInc(INT i)
{	SAFInc = i;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getSAFInc() const
{	return	SAFInc;	}

template <class TMOType>
void	ConfigurationSAFNr<TMOType>::setEnergy(RootEnergies e)
{	rootEnergies = e;	}

template <class TMOType>
const RootEnergies &	ConfigurationSAFNr<TMOType>::getEnergy() const
{	return	rootEnergies;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getNExt() const
{	return	nExtOpen + (nExtClosed << 1);	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getNExtOpen() const
{	return	nExtOpen;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getNExtClosed() const
{	return	nExtClosed;	}

template <class TMOType>
Configuration<TMOType>	ConfigurationSAFNr<TMOType>::getInternal() const
{	return	internal;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::isReference() const
{	return isRef;	}

template <class TMOType>
INT	ConfigurationSAFNr<TMOType>::getSelectionFlags() const
{	return SelectionFlags;	}


template <class TMOType>
void	ConfigurationSAFNr<TMOType>::setInternal(Configuration<TMOType> conf)
{	internal = conf;	}

//--------------------------------------------------------------------------


template <class TMOType>
ostream& operator<<(ostream& s, const ConfigurationSAFNr<TMOType> & conf)
{
	if ( conf.isReference() )
		s << "  ref.  ";
	else
	if ( conf.getSelectionFlags() & 1 )
		s << " degen. ";
	else
	if ( conf.getSelectionFlags() & 2 )
		s << " uncond.";
	else
		s << "        ";
	s << setprecision(10) << setiosflags(ios::fixed) << conf.getEnergy() << " " 
		<< *((Configuration<TMOType> *) &conf);

//	s <<  conf.getSAFNr() << " "
//		<< conf.getSAFInc() << endl;

	return s;
}


/*
istream& operator>>(istream & s, ConfigurationSAFNr<MOType> &conf)
{
	return s;
}
*/

//--------------------------------------------------------------------------

template class ConfigurationSAFNr<MOType>;
template class ConfigurationSAFNr<GeneralizedMO>;

template ostream& operator << (ostream& s, const
    ConfigurationSAFNr<GeneralizedMO> & v);
template ostream& operator << (ostream& s, const
    ConfigurationSAFNr<MOType> & v);


//**********************************************************************
//	explicit instantiation for g++ 2.7.2 
//	("template <class TMOType> void calcInteraction(...)" does not work)
/*
void	ConfigurationSAFNrInstantiateTemplates()
{
	{
	ConfigurationSAFNr<MOType> *a = NULL;
		cout << *a;
//		cin >> *a;
	}

	{
	Configuration<GeneralizedMO> *a = NULL;
		cout << *a;
	}
}
*/
