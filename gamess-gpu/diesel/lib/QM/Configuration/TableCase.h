//***********************************************************************
//
//	Name:	TableCases.h
//
//	Description:	implements configuration handling
//					based on second quantization
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	08.06.1996
//
//	Literature:	Volker Pleﬂ: "Ein direktes, individuell
//			selektierendes Multireferenz Konfigurations-
//			Wechselwirkungsverfahren", 
//			Dissertation Uni Bonn, 1994
//
//
//
//***********************************************************************

#ifndef __TABLECASE_H
#define __TABLECASE_H

#include "../../../config.h"

#include "Configuration.h"

#include <iostream>

#include "../../Math/etc/MathObject.h"
#include "../MO/GeneralizedMO.h"
#include "DiffConf.h"
#include "../../Math/etc/BinomialCoefficient.h"
#include "../IntegralIndex/TwoElectronIntegralIndex.h"
#include "../IntegralIndex/TwoElectronIntegralCbExIndex.h"

class	MOContainer;

template <class TMOType> class TableCase;
template <class TMOType> ostream & operator<< (ostream & s, const TableCase<TMOType> &);

template <class TMOType>
class TableCase : public MathObject {
public:
	TableCase(BinomialCoefficient *);
	TableCase(
		INT _openShells, INT Nr);
	TableCase(
		INT _openShells, INT _dk, INT __P, INT __R, INT _qR = 1, INT _qL = 1);
	TableCase(BinomialCoefficient *, 
		INT _openShells, INT _dk, INT __P, INT __R, INT _qR = 1, INT _qL = 1);
	~TableCase() {}
	
//----------------------------------------------------------------------	
	
	INT	calc(const Configuration<TMOType> & a,
		const Configuration<TMOType> & b);

	INT	calc(const DiffConf<TMOType> & dc);

	void	calcQ(const Configuration<TMOType> & a,
		const Configuration<TMOType> & b);

	void	calcQ(const DiffConf<TMOType> & dc);
	
	void	updateqR(const DiffConf<TMOType> & dc);
	
	
	
	// the following only for configurations which are known
	// to be less or equal an double excitation (faster code)
	void	calcLess3NoInd(const Configuration<TMOType> & a,
		const Configuration<TMOType> & b);

	void	calcLess3NoInd(const DiffConf<TMOType> & dc);


	// additionally calculate integral indices
	void	calcLess3(
		const DiffConf<TMOType> & dc);

	void	calcLess3(const DiffConf<TMOType> & dc,
		TwoElectronIntegralIndex<TMOType> & Cb,
		TwoElectronIntegralIndex<TMOType> & Ex);

		
//----------------------------------------------------------------------	
	
	
	void	setNr(INT _openShells, INT Nr, INT _qR = 1, INT _qL = 1);
	const char	*getName() const;
	INT	getNumberOfMoreOpenShells() const;
	INT	getdK() const;
	INT	getP() const;
	INT	getR() const;
	INT	getqR() const;
	INT	getqL() const;
	
	const TwoElectronIntegralCbExIndex & getCbExIndex() const;
	TwoElectronIntegralCbExIndex & getCbExIndex();

//----------------------------------------------------------------------	
	
	friend ostream & operator<< <TMOType> (ostream & s, const TableCase<TMOType> &);


	INT	operator == (const TableCase<TMOType> &) const;
	INT	operator != (const TableCase<TMOType> &) const;

	
private:
	INT calcQfromPos(const INT *ipos, INT npos);
	

BinomialCoefficient	*binom;
INT	openShells;	// number of open shells in conf. of higher super categorie
INT	dK;			// difference between supercategories ("delta K")
INT	P;			// type of integrals
INT	R;			// classification of interaction
INT	qR;			// classification of interacting open shells 
				// in "right" configuration
INT	qL;			// classification of interacting open shells
				// in "left" configuration
TwoElectronIntegralCbExIndex CbExIndex;
};

template <class TMOType>	ostream & operator<<(ostream & s, const TableCase<TMOType> &);


//=========================================================================


template <class TMOType>
inline
TableCase<TMOType>::TableCase(BinomialCoefficient *_binom)
{
	binom = _binom;
	openShells = dK = P = R = 0;
	qR = 0;
	qL = 0;
}

template <class TMOType>
inline
INT TableCase<TMOType>::calcQfromPos(const INT *ipos, INT npos)
{
INT	q = 1;

//	cout << "npos= " << npos << endl;

	for ( INT i=0 ; i<npos ; i++ )
	{//	cout << "i, *ipos = " << i << ", " << *ipos << endl;
		q += binom->get(*ipos++, i+1);
	}
//	cout << "q= " << q << endl;	

	return q;
}


template <class TMOType>
inline
TableCase<TMOType>::TableCase(BinomialCoefficient *_binom,
	INT _openShells, INT _dK, INT __P, INT __R, INT _qR, INT _qL)
{
	binom = _binom;
	openShells = _openShells;
	dK = _dK;
	P = __P;
	R = __R;
	qR = _qR;
	qL = _qL;
	calcQfromPos((const INT *) &dK, 0);
}

template <class TMOType>
inline
TableCase<TMOType>::TableCase(
	INT _openShells, INT _dK, INT __P, INT __R, INT _qR, INT _qL)
{
	openShells = _openShells;
	dK = _dK;
	P = __P;
	R = __R;
	qR = _qR;
	qL = _qL;
}

template <class TMOType>
inline
INT	TableCase<TMOType>::getNumberOfMoreOpenShells() const
{	return openShells;	}

template <class TMOType>
inline
INT	TableCase<TMOType>::getdK() const
{	return dK;	}

template <class TMOType>
inline
INT	TableCase<TMOType>::getP() const
{	return P;	}

template <class TMOType>
inline
INT	TableCase<TMOType>::getR() const
{	return R;	}

template <class TMOType>
inline
INT	TableCase<TMOType>::getqR() const
{	return qR;	}

template <class TMOType>
inline
INT	TableCase<TMOType>::getqL() const
{	return qL;	}


template <class MOType>
inline
const TwoElectronIntegralCbExIndex & TableCase<MOType>::getCbExIndex() const
{	return	CbExIndex;	}

template <class MOType>
inline
TwoElectronIntegralCbExIndex & TableCase<MOType>::getCbExIndex()
{	return	CbExIndex;	}

template <class TMOType>
inline
void	TableCase<TMOType>::updateqR(const DiffConf<TMOType> & dc)
{
	if ( dc.getTo().getNumberOfOpenShells() )
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());
}

template <class TMOType>
inline
INT	TableCase<TMOType>::operator == (const TableCase<TMOType> & t) const
{
	if ( dK!=t.dK )
		return 0;	
	if ( R!=t.R )
		return 0;	
	if ( P!=t.P )
		return 0;	
	if ( qR!=t.qR )
		return 0;	
	if ( qL!=t.qL )
		return 0;	
	return 1;
}

template <class TMOType>
inline
INT	TableCase<TMOType>::operator != (const TableCase<TMOType> & t) const
{
	if ( dK!=t.dK )
		return 1;	
	if ( R!=t.R )
		return 1;	
	if ( P!=t.P )
		return 1;	
	if ( qR!=t.qR )
		return 1;	
	if ( qL!=t.qL )
		return 1;	
	return 0;
}


#endif
