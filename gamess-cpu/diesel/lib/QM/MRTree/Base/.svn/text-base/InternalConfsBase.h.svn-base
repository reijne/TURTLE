//***********************************************************************
//
//	Name:			InternalConfsBase.h
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

#ifndef __InternalConfsBase_H
#define __InternalConfsBase_H

#include "../../../../config.h"

#include <iostream>
using std::ostream;
using std::istream;

class	NExternalsSet;
class	TupelStructureSet;

class InternalConfsBase {
public:
	InternalConfsBase(INT nExt);
	InternalConfsBase(istream &s);


	INT	getNExt() const;

	//----------------------------------------------------------------	

	INT operator == (const InternalConfsBase &) const;
	INT operator != (const InternalConfsBase &) const;
	INT operator <= (const InternalConfsBase &) const;

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const InternalConfsBase &);
	friend istream& operator>>(istream & s, InternalConfsBase &);

	//----------------------------------------------------------------	
protected:
INT	nExt;					// number externals
};


inline
INT	InternalConfsBase::getNExt() const
{	return nExt;	}


inline
INT InternalConfsBase::operator == (const InternalConfsBase &a) const
{	return nExt==a.getNExt();	}

inline
INT InternalConfsBase::operator != (const InternalConfsBase &a) const
{	return nExt!=a.getNExt();	}

inline
INT InternalConfsBase::operator <= (const InternalConfsBase &a) const
{	return nExt<=a.getNExt();	}


#endif
