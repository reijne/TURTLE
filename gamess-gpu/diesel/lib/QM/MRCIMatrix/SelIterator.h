//***********************************************************************
//
//	Name:			SelIterator.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			28.10.1998
//
//
//
//
//
//***********************************************************************


#ifndef __SelIterator_h
#define __SelIterator_h

#include "../../../config.h"


#include "../MRTree/Container/ContainerIterator.h"

#include "../MO/MOType.h"

template <class T> class AVLSet;
template <class T> class ConfigurationStart;
class NExternalsDiag;
class InternalConfsDiag;
class TupelStructureDiag;
class extMOsDiag;
class MOAccess;

class SelIterator {
public:
	SelIterator(
		const NExternalsDiag & Psi,
		const AVLSet<ConfigurationStart<MOType> > *internal,
		MOAccess	* const *,
		INT maxExc
		);
	~SelIterator();

	INT	isEnd() const;
	void	next(INT first = 0);
	
	INT	getConfStart() const;
	INT getSAFStart() const;
	INT	getSAFInc() const;
	

private:
	void	prepair(INT);


const NExternalsDiag &Psi;
const AVLSet<ConfigurationStart<MOType> > *internal;
MOAccess	* const *moAccess;
INT	maxExc;


INT nExt;
ContainerIterator	iter0;
ContainerIterator	iter1;
ContainerIterator	iter2;
INT	iter3;
const NExternalsDiag *p0;
const InternalConfsDiag	*p1;
const TupelStructureDiag	*p2;
const extMOsDiag	*p3;
const MOType	*mo;

INT	lconfStart;
INT	lSAFStart;
INT	internalOpen;

INT	end;
INT	confStart;
INT	SAFStart;
INT	SAFInc;
};





inline
INT	SelIterator::isEnd() const
{	return	end;	}

inline
INT	SelIterator::getConfStart() const
{	return confStart;	}

inline
INT SelIterator::getSAFStart() const
{	return SAFStart;	}

inline
INT	SelIterator::getSAFInc() const
{	return SAFInc;	}


#endif



