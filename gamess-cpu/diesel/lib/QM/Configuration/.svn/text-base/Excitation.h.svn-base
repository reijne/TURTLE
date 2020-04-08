//***********************************************************************
//
//	Name:	Excitation.h
//
//	Description:	general representation of excitations
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	11.06.1996
//
//
//
//
//
//***********************************************************************

#ifndef __EXCITATION_H
#define __EXCITATION_H

#include "../../../config.h"

#include <iostream>
#include <stdarg.h>

#include "../../Math/etc/MathObject.h"

#include "ConfigurationGlobals.h"
#include "../MO/MOType.h"

template <class TMOType> class	Configuration;


class Excitation : public MathObject {
public:
	Excitation();
	Excitation(INT Order, MOType from ...);
	~Excitation() {}
	
	
	// Copy-Konstruktor
	Excitation(const Excitation &);
	
	// Zuweisungs-Operator
	Excitation & operator = (const Excitation &);


//--------------------------------------------------------------------------

	void	calcExcitation(const Configuration<MOType> & a,
		const Configuration<MOType> & b);
//	void	calcExcitation(const Configuration & a,
//		const Configuration & b);

//--------------------------------------------------------------------------


	friend INT operator == (const Excitation &, const Excitation &);
	friend INT operator != (const Excitation &, const Excitation &);

//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const Excitation &);

//--------------------------------------------------------------------------

	INT	getOrder() const;
	void	setOrder(INT n);

	const MOType	*getFromP() const;
	MOType	getFrom(INT n) const;
	void	setFrom(INT n, MOType i);

	const MOType	*getToP() const;
	MOType	getTo(INT n) const;
	void	setTo(INT n, MOType i);
	

	INT	operator < (const Excitation &exc) const
	{
		if ( order<exc.order )
			return 1;
		if ( order>exc.order )
			return 0;
		for ( INT i=0 ; i<order ; i++ )
		{
			if ( pFrom[i]<exc.pFrom[i] )
				return 1;
			if ( pFrom[i]>exc.pFrom[i] )
				return 0;
			if ( pTo[i]<exc.pTo[i] )
				return 1;
			if ( pTo[i]>exc.pTo[i] )
				return 0;
		}
		return 0;
	}

private:
INT	order;
MOType	pFrom[MAXELECTRONS];
MOType	pTo[MAXELECTRONS];
};



inline
Excitation::Excitation()
{	order = 0;	}


inline
INT	Excitation::getOrder() const
{	return	order;	}

inline
void	Excitation::setOrder(INT n)
{	order = n;	}

inline
const MOType	*Excitation::getFromP() const
{	return	pFrom;	}

inline
MOType	Excitation::getFrom(INT n) const
{	return	pFrom[n];	}

inline
void	Excitation::setFrom(INT n, MOType i)
{	pFrom[n] = i;	}

inline
const MOType	*Excitation::getToP() const
{	return	pTo;	}

inline
MOType	Excitation::getTo(INT n) const
{	return	pTo[n];	}

inline
void	Excitation::setTo(INT n, MOType i)
{	pTo[n] = i;	}


#endif
